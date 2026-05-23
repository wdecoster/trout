use dbscan::{Classification, Model};

/// Returns a boolean mask: true = noise (outlier), false = in a cluster.
pub fn find_outliers(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<bool> {
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
