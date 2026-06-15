use dbscan::{Classification, Model};

/// Returns a boolean mask: true = noise (outlier), false = in a cluster.
///
/// Implementation: custom DBSCAN over a precomputed pairwise squared-distance matrix. The
/// `dbscan` crate version is kept under [`find_outliers_crate`] for benchmarking and as a
/// reference oracle; the comparison tests below assert that the two agree on noise/non-noise
/// classification (the cluster IDs themselves can legitimately differ between implementations
/// because of border-point assignment ambiguity).
#[cfg(test)]
pub fn find_outliers(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<bool> {
    find_clusters(points, eps, min_samples)
        .iter()
        .map(|&label| label < 0)
        .collect()
}

/// DBSCAN cluster labels, one per point: `-1` = noise (outlier), `>= 0` = cluster id. Same
/// clustering as [`find_outliers`], but keeps the cluster identity so callers can, e.g., pick the
/// cluster with the longest alleles as a reference. Cheap: the cluster expansion already walks
/// each connected component, so retaining its id costs nothing extra.
pub fn find_clusters(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<i32> {
    if points.len() < min_samples {
        // Too few points to form a cluster; match find_outliers' "nobody is an outlier" by
        // placing everyone in one trivial cluster.
        return vec![0; points.len()];
    }
    find_outliers_matrix(points, eps, min_samples)
}

/// Custom DBSCAN: precomputes the n×n distance matrix once (upper triangle only), then runs
/// standard cluster expansion against it. Wins over the brute-force `dbscan` crate by:
///   (a) computing each pairwise distance exactly once instead of once per range query, and
///   (b) avoiding per-range-query Vec allocations (the flamegraph showed those dominated).
fn find_outliers_matrix(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<i32> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }

    let eps_sq = eps * eps;

    // Collect all neighbour edges in a single flat Vec first, then build CSR-style adjacency
    // from a single allocation pair. This removes the ~12% Vec-realloc cost the flamegraph
    // showed coming from growing 1000 separate small Vec<usize>'s with capacity doubling.
    // Self-inclusion: standard DBSCAN counts the point itself in its own ε-neighbourhood, so
    // a "core" point has at least min_samples members including itself — we add self-loops
    // when building the adjacency below.
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for i in 0..n {
        let pi = points[i].as_slice();
        for (j, pj_vec) in points.iter().enumerate().skip(i + 1) {
            let pj = pj_vec.as_slice();
            // Iterator-zipped sum: LLVM auto-vectorises this on f64 slices and elides bounds
            // checks. This is the inner-most hot loop — it runs ~n²/2 × dim times per locus.
            let dsq: f64 = pi
                .iter()
                .zip(pj.iter())
                .map(|(a, b)| {
                    let d = a - b;
                    d * d
                })
                .sum();
            if dsq <= eps_sq {
                edges.push((i as u32, j as u32));
            }
        }
    }

    // Build CSR adjacency: row_starts[i..i+1] indexes col_indices[] for point i's neighbours.
    let mut row_lens = vec![1usize; n]; // each point includes itself
    for &(i, j) in &edges {
        row_lens[i as usize] += 1;
        row_lens[j as usize] += 1;
    }
    let mut row_starts = vec![0usize; n + 1];
    for i in 0..n {
        row_starts[i + 1] = row_starts[i] + row_lens[i];
    }
    let total = row_starts[n];
    let mut col_indices: Vec<u32> = vec![0u32; total];
    // Cursors track the next free slot per row as we fill.
    let mut cursor = row_starts.clone();
    for i in 0..n {
        col_indices[cursor[i]] = i as u32; // self-loop
        cursor[i] += 1;
    }
    for &(i, j) in &edges {
        col_indices[cursor[i as usize]] = j;
        cursor[i as usize] += 1;
        col_indices[cursor[j as usize]] = i;
        cursor[j as usize] += 1;
    }
    drop(edges);
    drop(cursor);
    drop(row_lens);

    // Helper to read row i's neighbours as a slice.
    let neighbours_of = |i: usize| -> &[u32] { &col_indices[row_starts[i]..row_starts[i + 1]] };

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum St {
        Unvisited,
        Noise,
        Clustered,
    }
    let mut status = vec![St::Unvisited; n];
    // Cluster id per point; -1 = noise/unassigned. Filled alongside `status` as the BFS walks
    // each connected component, so callers can identify individual clusters.
    let mut labels = vec![-1i32; n];
    let mut cluster_id = 0i32;
    let mut queue: Vec<usize> = Vec::new();

    for p in 0..n {
        if status[p] != St::Unvisited {
            continue;
        }
        if neighbours_of(p).len() < min_samples {
            status[p] = St::Noise;
            continue;
        }
        // p is a core point — start a new cluster and seed the expansion frontier with
        // its neighbours. We mark a point Clustered *when we enqueue it*, not when we pop,
        // so duplicate enqueues are impossible without a separate "in-queue" mask.
        status[p] = St::Clustered;
        labels[p] = cluster_id;
        queue.clear();
        for &q in neighbours_of(p) {
            let q = q as usize;
            if q != p && (status[q] == St::Unvisited || status[q] == St::Noise) {
                status[q] = St::Clustered;
                labels[q] = cluster_id;
                queue.push(q);
            }
        }

        let mut idx = 0;
        while idx < queue.len() {
            let q = queue[idx];
            idx += 1;
            // Only core points propagate the cluster (border points are absorbed but don't
            // expand further).
            if neighbours_of(q).len() >= min_samples {
                for &r in neighbours_of(q) {
                    let r = r as usize;
                    if status[r] == St::Unvisited || status[r] == St::Noise {
                        status[r] = St::Clustered;
                        labels[r] = cluster_id;
                        queue.push(r);
                    }
                }
            }
        }
        cluster_id += 1;
    }

    labels
}

/// Reference implementation using the `dbscan` crate. Kept for comparison tests and for any
/// future bench harness. Not on any hot path.
#[allow(dead_code)]
fn find_outliers_crate(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<bool> {
    if points.len() < min_samples {
        return vec![false; points.len()];
    }
    let model = Model::new(eps, min_samples);
    let population: Vec<Vec<f64>> = points.to_vec();
    model
        .run(&population)
        .into_iter()
        .map(|c| matches!(c, Classification::Noise))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_same(points: &[Vec<f64>], eps: f64, ms: usize) {
        // find_outliers_matrix now returns cluster labels; compare noise classification (label < 0).
        let custom: Vec<bool> = find_outliers_matrix(points, eps, ms)
            .iter()
            .map(|&label| label < 0)
            .collect();
        let crate_v = find_outliers_crate(points, eps, ms);
        assert_eq!(
            custom, crate_v,
            "noise classification disagreement at eps={}, min_samples={}",
            eps, ms
        );
    }

    fn build_clusters(seed: u64) -> Vec<Vec<f64>> {
        // Deterministic 3 dense clusters in 2D + scattered noise. The seed permutes
        // tie-breaker order so border-point classification varies between calls; the noise
        // mask should be stable regardless.
        let mut points: Vec<Vec<f64>> = Vec::new();
        let centers = [(0.0, 0.0), (10.0, 0.0), (5.0, 8.0)];
        for &(cx, cy) in &centers {
            for i in 0..25 {
                let a = (i as f64 + seed as f64) * 0.37;
                points.push(vec![cx + a.sin() * 0.5, cy + a.cos() * 0.5]);
            }
        }
        // Eight scattered noise points well outside any eps neighbourhood.
        for i in 0..8 {
            let s = i as f64;
            points.push(vec![-50.0 + s * 17.0, 100.0 - s * 13.0]);
        }
        points
    }

    #[test]
    fn equivalence_simple_clusters() {
        let pts = build_clusters(0);
        assert_same(&pts, 1.0, 5);
    }

    #[test]
    fn equivalence_varied_params() {
        let pts = build_clusters(1);
        for eps in [0.3, 0.7, 1.5, 3.0] {
            for ms in [3, 5, 10] {
                assert_same(&pts, eps, ms);
            }
        }
    }

    #[test]
    fn equivalence_all_noise() {
        // Points so spread that nothing clusters at this eps.
        let pts: Vec<Vec<f64>> = (0..20).map(|i| vec![(i * 100) as f64, 0.0]).collect();
        assert_same(&pts, 1.0, 5);
    }

    #[test]
    fn equivalence_one_big_cluster() {
        // Everything within eps of everything else.
        let pts: Vec<Vec<f64>> = (0..30).map(|i| vec![i as f64 * 0.01, 0.0]).collect();
        assert_same(&pts, 1.0, 5);
    }

    #[test]
    fn equivalence_high_dim() {
        // Stress at k=5-style dimensionality (~209 dims) on a small point set.
        let dim = 200;
        let mut pts: Vec<Vec<f64>> = Vec::new();
        for c in 0..3 {
            for i in 0..15 {
                let mut v = vec![0.0; dim];
                v[c] = 1.0 + (i as f64 * 0.001);
                pts.push(v);
            }
        }
        assert_same(&pts, 0.5, 5);
    }

    #[test]
    fn empty_input() {
        let pts: Vec<Vec<f64>> = Vec::new();
        let r = find_outliers(&pts, 0.5, 5);
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn fewer_points_than_min_samples_returns_all_false() {
        // Trout treats "too few points to analyse" as "no outliers", not "all outliers".
        let pts: Vec<Vec<f64>> = (0..3).map(|i| vec![i as f64, 0.0]).collect();
        let r = find_outliers(&pts, 0.1, 5);
        assert_eq!(r, vec![false; 3]);
    }
}
