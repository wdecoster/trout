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

pub fn render_scatter_plots(data: &[LocusPlotData], path: &Path) {
    if data.is_empty() {
        eprintln!("No outlier loci to plot");
        return;
    }

    let n_data = data.len();
    // +1 for the shared legend panel
    let n_total = n_data + 1;
    let cols = n_total.min(4);
    let rows = (n_total + cols - 1) / cols;

    let (mut all_plots, mut all_layouts): (Vec<Vec<Plot>>, Vec<Layout>) = data
        .iter()
        .map(|locus| {
            let mut panel: Vec<Plot> = Vec::new();

            if !locus.normal.is_empty() {
                let xy: Vec<(f64, f64)> = locus.normal.iter().map(|(l, f, _)| (*l, *f)).collect();
                let labels: Vec<String> = locus.normal.iter().map(|(_, _, s)| s.clone()).collect();
                panel.push(Plot::Scatter(
                    ScatterPlot::new()
                        .with_data(xy)
                        .with_tooltip_labels(labels)
                        .with_color("steelblue")
                        .with_size(4.0)
                        .with_marker_opacity(0.6),
                ));
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
                        .with_color("crimson")
                        .with_size(6.0)
                        .with_marker_opacity(0.8),
                ));
            }

            let mut layout = Layout::auto_from_plots(&panel)
                .with_title(&locus.title)
                .with_x_label("length (bp)")
                .with_y_label(&locus.y_label)
                .with_interactive();

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
    let legend_panel = vec![
        Plot::Scatter(
            ScatterPlot::new()
                .with_data(vec![(1.0, 3.0)])
                .with_color("steelblue")
                .with_size(6.0)
                .with_marker_opacity(0.6)
                .with_legend("normal"),
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

    let scene = Figure::new(rows, cols)
        .with_plots(all_plots)
        .with_layouts(all_layouts)
        .render();

    let svg = SvgBackend::new().render_scene(&scene);
    std::fs::write(path, svg).expect("Failed to write plot SVG");
    eprintln!("Wrote scatter plot ({} loci) to {}", n_data, path.display());
}
