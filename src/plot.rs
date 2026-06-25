use kuva::backend::svg::SvgBackend;
use kuva::plot::scatter::ScatterPlot;
use kuva::render::annotations::TextAnnotation;
use kuva::render::figure::Figure;
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use std::path::Path;

pub struct LocusPlotData {
    pub title: String,
    pub y_label: String,
    pub normal: Vec<(f64, f64, String)>, // (length_bp, kmer_freq, sample)
    pub outliers: Vec<(f64, f64, String)>, // outliers for this axis
    pub other_outliers: Vec<(f64, f64, String)>, // outliers flagged on a different axis
}

/// Collapse points that coincide at plotting resolution into unique positions,
/// counting how many samples fell on each. The normal cloud has one entry per
/// sample but most overlap heavily; deduping to unique (length-to-nearest-bp,
/// frequency-to-1e-3) positions cuts the marker count to the number of distinct
/// values, while the count lets us shade dense spots so a position holding 300
/// samples no longer looks like a single one. Returns `(x, y, count)`.
fn dedup_with_counts(points: &[(f64, f64, String)]) -> Vec<(f64, f64, usize)> {
    use std::collections::HashMap;
    let mut idx: HashMap<(i64, i64), usize> = HashMap::new();
    let mut out: Vec<(f64, f64, usize)> = Vec::new();
    for (x, y, _) in points {
        let key = (x.round() as i64, (y * 1000.0).round() as i64);
        match idx.get(&key) {
            Some(&i) => out[i].2 += 1,
            None => {
                idx.insert(key, out.len());
                out.push((*x, *y, 1));
            }
        }
    }
    out
}

/// Map a sample count to a blue shade — light for few, dark for many — so the
/// density a single deduped marker hides is still visible. Log-scaled because the
/// counts are heavily skewed; `max_count` is the densest position in the same
/// panel, so shading is relative within each plot.
fn count_color(count: usize, max_count: usize) -> String {
    let t = if max_count <= 1 {
        0.0
    } else {
        (count as f64).ln() / (max_count as f64).ln()
    };
    // Interpolate light steelblue -> dark navy.
    let lerp = |a: f64, b: f64| (a + (b - a) * t).round() as i64;
    format!(
        "rgb({}, {}, {})",
        lerp(158.0, 8.0),
        lerp(202.0, 48.0),
        lerp(225.0, 107.0)
    )
}

pub fn render_scatter_plots(data: &[LocusPlotData], path: &Path, jitter: bool, interactive: bool) {
    if data.is_empty() {
        // Still write the output file so workflows (e.g. Snakemake) that declare the
        // --plot path as an expected output don't fail when there are no outliers.
        write_placeholder_svg(path);
        return;
    }

    // Only the scalable default collapses the normal cloud. Both --interactive
    // (every point must stay individually hoverable) and --jitter (points are
    // fanned out to show density by spread) need every point, so they skip dedup.
    let dedup = !interactive && !jitter;

    let n_data = data.len();
    // +1 for the shared legend panel
    let n_total = n_data + 1;
    let cols = n_total.min(4);
    let rows = n_total.div_ceil(cols);

    let (mut all_plots, mut all_layouts): (Vec<Vec<Plot>>, Vec<Layout>) = data
        .iter()
        .map(|locus| {
            let mut panel: Vec<Plot> = Vec::new();

            if !locus.normal.is_empty() {
                // The normal "cloud" is the bulk of the points and is just context.
                let normal_scatter = if dedup {
                    // Default: collapse coincident points (most of the cloud
                    // overlaps) so the series renders as batched plain <circle>s
                    // instead of one <g><title><circle> per sample, and shade each
                    // kept position by how many samples it holds so the density
                    // survives the dedup.
                    let counted = dedup_with_counts(&locus.normal);
                    let max_c = counted.iter().map(|(_, _, c)| *c).max().unwrap_or(1);
                    let xy: Vec<(f64, f64)> = counted.iter().map(|(x, y, _)| (*x, *y)).collect();
                    let colors: Vec<String> = counted
                        .iter()
                        .map(|(_, _, c)| count_color(*c, max_c))
                        .collect();
                    ScatterPlot::new()
                        .with_data(xy)
                        .with_colors(colors)
                        .with_size(4.0)
                        .with_marker_opacity(0.85)
                } else {
                    // --interactive or --jitter: keep every point. Under
                    // --interactive, attach the sample name so any point is
                    // hoverable; otherwise overlap/opacity conveys local density.
                    let xy: Vec<(f64, f64)> =
                        locus.normal.iter().map(|(x, y, _)| (*x, *y)).collect();
                    let s = ScatterPlot::new()
                        .with_data(xy)
                        .with_color("steelblue")
                        .with_size(4.0)
                        .with_marker_opacity(0.6);
                    if interactive {
                        let labels: Vec<String> =
                            locus.normal.iter().map(|(_, _, n)| n.clone()).collect();
                        s.with_tooltip_labels(labels)
                    } else {
                        s
                    }
                };
                panel.push(Plot::Scatter(normal_scatter));
            }

            if !locus.other_outliers.is_empty() {
                let xy: Vec<(f64, f64)> = locus
                    .other_outliers
                    .iter()
                    .map(|(l, f, _)| (*l, *f))
                    .collect();
                let labels: Vec<String> = locus
                    .other_outliers
                    .iter()
                    .map(|(_, _, s)| s.clone())
                    .collect();
                panel.push(Plot::Scatter(
                    ScatterPlot::new()
                        .with_data(xy)
                        .with_tooltip_labels(labels)
                        .with_tooltips()
                        .with_color("darkorange")
                        .with_size(5.0)
                        .with_marker_opacity(0.8),
                ));
            }

            if !locus.outliers.is_empty() {
                let xy: Vec<(f64, f64)> = locus.outliers.iter().map(|(l, f, _)| (*l, *f)).collect();
                let labels: Vec<String> =
                    locus.outliers.iter().map(|(_, _, s)| s.clone()).collect();
                panel.push(Plot::Scatter(
                    ScatterPlot::new()
                        .with_data(xy)
                        .with_tooltip_labels(labels)
                        .with_tooltips()
                        .with_color("crimson")
                        .with_size(6.0)
                        .with_marker_opacity(0.8),
                ));
            }

            // By default NOT interactive: kuva's interactive mode wraps every point
            // in a <g class="tt"><title> group and disables the batched CircleBatch
            // fast path, which is what made these SVGs explode. Outlier names are
            // shown via the arrow annotations below, and a native <title> (from
            // `with_tooltips()` on the outlier series) still gives hover. With
            // --interactive the user has opted into the heavy per-point machinery.
            let mut layout = Layout::auto_from_plots(&panel)
                .with_title(&locus.title)
                .with_x_label("length (bp)")
                .with_y_label(&locus.y_label);
            if interactive {
                layout = layout.with_interactive();
            }

            // Place each outlier label inside the data range.
            // Direction is based on which side of the data midpoint the point sits on —
            // unlike centroid-based shifting, this is invariant to skewed distributions
            // where a dense cluster pulls the centroid past the outlier point.
            let all_x: Vec<f64> = locus
                .normal
                .iter()
                .chain(locus.outliers.iter())
                .chain(locus.other_outliers.iter())
                .map(|(x, _, _)| *x)
                .collect();
            let all_y: Vec<f64> = locus
                .normal
                .iter()
                .chain(locus.outliers.iter())
                .chain(locus.other_outliers.iter())
                .map(|(_, y, _)| *y)
                .collect();
            let x_min = all_x.iter().cloned().fold(f64::INFINITY, f64::min);
            let x_max = all_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let y_min = all_y.iter().cloned().fold(f64::INFINITY, f64::min);
            let y_max = all_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let x_range = x_max - x_min;
            let y_range = y_max - y_min;
            let x_mid = (x_min + x_max) / 2.0;
            let y_mid = (y_min + y_max) / 2.0;
            let off_frac = 0.15;

            // Compute initial label positions, then stagger any that overlap in y.
            let mut labels: Vec<(f64, f64, f64, f64, String)> = locus
                .outliers
                .iter()
                .map(|(px, py, sample)| {
                    let lx = if *px <= x_mid {
                        px + off_frac * (x_range + 1e-10)
                    } else {
                        px - off_frac * (x_range + 1e-10)
                    };
                    let ly = if *py <= y_mid {
                        py + off_frac * (y_range + 1e-10)
                    } else {
                        py - off_frac * (y_range + 1e-10)
                    };
                    (*px, *py, lx, ly, sample.clone())
                })
                .collect();

            // Sort by label position, then push overlapping labels apart in y.
            labels.sort_by(|a, b| {
                a.2.partial_cmp(&b.2)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
            });
            let y_step = 0.07 * (y_range + 1e-10);
            for i in 1..labels.len() {
                let prev_ly = labels[i - 1].3;
                if (labels[i].3 - prev_ly).abs() < y_step {
                    labels[i].3 = prev_ly + y_step;
                }
            }

            for (px, py, lx, ly, sample) in &labels {
                layout = layout.with_annotation(
                    TextAnnotation::new(sample, *lx, *ly)
                        .with_arrow(*px, *py)
                        .with_font_size(9),
                );
            }

            (panel, layout)
        })
        .unzip();

    // Shared legend panel: one dot per series category so the legend appears once.
    // The density-shading hint only applies in the default (deduped) mode.
    let normal_legend = if dedup {
        "normal (darker = more samples)"
    } else {
        "normal"
    };
    let legend_panel = vec![
        Plot::Scatter(
            ScatterPlot::new()
                .with_data(vec![(1.0, 3.0)])
                .with_color("steelblue")
                .with_size(6.0)
                .with_marker_opacity(0.85)
                .with_legend(normal_legend),
        ),
        Plot::Scatter(
            ScatterPlot::new()
                .with_data(vec![(1.0, 2.0)])
                .with_color("darkorange")
                .with_size(7.0)
                .with_marker_opacity(0.8)
                .with_legend("outlier (other axis)"),
        ),
        Plot::Scatter(
            ScatterPlot::new()
                .with_data(vec![(1.0, 1.0)])
                .with_color("crimson")
                .with_size(8.0)
                .with_marker_opacity(0.8)
                .with_legend("outlier"),
        ),
    ];
    let legend_layout = Layout::auto_from_plots(&legend_panel)
        .with_title("Legend")
        .with_x_label("")
        .with_y_label("");
    all_plots.push(legend_panel);
    all_layouts.push(legend_layout);

    let mut scene = Figure::new(rows, cols)
        .with_plots(all_plots)
        .with_layouts(all_layouts)
        .render();
    // Figure::render only propagates per-panel interactivity (the point tooltips),
    // not the scene-level flag that emits the search/save UI strip — so set it
    // explicitly when the user asked for full interactivity.
    if interactive {
        scene.interactive = true;
    }

    let svg = SvgBackend::new().render_scene(&scene);
    std::fs::write(path, svg).expect("Failed to write plot SVG");
    eprintln!("Wrote scatter plot ({} loci) to {}", n_data, path.display());
}

/// Write a minimal valid SVG carrying a "no outliers" message. Used when there are no
/// outlier loci to plot, so the --plot output file still exists for downstream workflows.
fn write_placeholder_svg(path: &Path) {
    let svg = r##"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="480" height="160" viewBox="0 0 480 160">
  <rect width="480" height="160" fill="white"/>
  <text x="240" y="80" font-family="sans-serif" font-size="18" fill="#555555" text-anchor="middle" dominant-baseline="middle">No outlier loci detected</text>
</svg>
"##;
    std::fs::write(path, svg).expect("Failed to write placeholder plot SVG");
    eprintln!(
        "No outlier loci to plot; wrote placeholder to {}",
        path.display()
    );
}
