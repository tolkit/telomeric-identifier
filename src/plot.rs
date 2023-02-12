use anyhow::Result;
use csv::ReaderBuilder;
use serde::Deserialize;
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

/// The plot margins
const MARGIN: i32 = 40;

/// The entry point for `tidk plot`.
pub fn plot(matches: &clap::ArgMatches) -> Result<()> {
    // parse the command line options
    let tsv = matches.get_one::<PathBuf>("tsv").expect("errored by clap");
    // a bug here for manual input of chromosome cut-off which I can't figure out right now.
    let chromosome_cutoff = 0;
    let height_subplot = matches.get_one::<i32>("height").expect("defualted by clap");
    let width = matches.get_one::<i32>("width").expect("defaulted by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("defaulted by clap");

    // parse the tsv
    let parsed_tsv = parse_tsv(tsv.to_path_buf())?;

    // calculate the number of chromosomes to plot with the length cutoff
    let chromosome_number = chromosome_number(&parsed_tsv, chromosome_cutoff);

    // height of plot
    let height: i32 = height_subplot * chromosome_number as i32 + (2 * MARGIN);

    // generate the plot data (see struct PlotData)
    let plot_data = generate_plot_data(parsed_tsv, height, *width, *height_subplot);

    // filter the data based on the cutoff
    let plot_data_filtered: Vec<PlotData> = plot_data
        .into_iter()
        .filter(|x| x.max > chromosome_cutoff as usize)
        .collect();

    // make the writable svg file
    let out_filename = format!("{}.svg", output.display());
    let mut svg_file = File::create(out_filename)?;

    // construct the svg
    let svg = format!(
            "<?xml version='1.0' encoding='UTF-8'  standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'> \
                \
                 <style type='text/css'> \
                 .chromosome_line:hover {{ stroke-opacity: 1.0; stroke: crimson; stroke-width: 2; }} \
                 </style> \
                \
                 {} \
                 </svg>",
            width,
            height,
            add_all_path_elements(plot_data_filtered, *height_subplot as isize, *width)
        );

    svg_file.write_all(svg.as_bytes())?;

    Ok(())
}

/// Deserialise the TSV records into a struct.
#[derive(Debug, Deserialize)]
pub struct TelomericRepeatRecord {
    pub id: String,
    pub window: i32,
    pub forward_repeat_number: i32,
    pub reverse_repeat_number: i32,
    pub telomeric_repeat: String,
}

/// This deserializes a TSV to a [`Vec<TelomericRepeatRecord>`].
fn parse_tsv(path: PathBuf) -> Result<Vec<TelomericRepeatRecord>> {
    let mut tsv_reader = ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let mut plot_coords_vec = Vec::new();

    for result in tsv_reader.deserialize() {
        let record: TelomericRepeatRecord = result?;
        plot_coords_vec.push(record);
    }

    Ok(plot_coords_vec)
}

/// Takes the parsed TSV and the chromosome cutoff,
/// loops through file to find the lengths of all the
/// chromosomes (to the nearest window) and reports the
/// number of elements.
fn chromosome_number(parsed_tsv: &Vec<TelomericRepeatRecord>, chromosome_cutoff: i32) -> usize {
    // so we can break the loop
    let file_length = parsed_tsv.len();
    // the iteration of the loop
    let mut it = 0usize;
    // a mutable vector to calculate svg path attribute
    let mut max_sizes = Vec::new();

    loop {
        if it == file_length - 1 {
            break;
        }

        if parsed_tsv[it].id == parsed_tsv[it + 1].id {
            it += 1;
            continue;
        } else {
            if parsed_tsv[it].window > chromosome_cutoff {
                max_sizes.push(parsed_tsv[it].window);
            }
            it += 1;
        }
    }
    // add 1 as we can never compare the last entry.
    max_sizes.len() + 1
}

/// Scale a range [min, max] to custom range [a, b].
/// Our range will be [0, height of subplots].
fn scale_y(y: f64, a: f64, b: f64, min: f64, max: f64) -> f64 {
    (((b - a) * (y - min)) / (max - min)) + a
}

/// Make the SVG path elements:
/// `Mx,yLx1,y1Lx2,y2`...
///
/// This is also where the paths are scaled to
/// the plot width and subplot height.
fn make_path_element(
    path_vec: Vec<(i32, i32)>,
    x_max: usize,
    y_max: usize,
    height: i32,
    width: i32,
    height_per_plot: i32,
) -> Option<String> {
    let subplot_gap = 25.0;
    // need this here...
    if path_vec.is_empty() {
        return None;
    }

    // somehow going to have to scale the lines to fit on all graphs
    let mut path = String::new();
    // margin at either side of the plot
    let width_incl_margins = width - (MARGIN * 2);

    let x_bin: f64 = width_incl_margins as f64 / x_max as f64;

    // add the first move to point
    path += &format!(
        "M{},{}",
        MARGIN as f64 / 2f64, // + x_bin, don't think I need x_bin here?
        // y is height - repeat number, a = 0, b = height_per_plot
        // min is again zero (no negative repeats), max is greatest repeats per chromosome
        height as f64
            - MARGIN as f64
            - scale_y(
                path_vec[0].1 as f64,
                0.0,
                height_per_plot as f64,
                0.0,
                y_max as f64 + subplot_gap
            )
    );

    let mut bin: f64 = x_bin;
    for element in path_vec.iter().skip(1) {
        path += &format!(
            "L{},{}",
            bin + (MARGIN as f64 / 2.0),
            //
            height as f64
                - MARGIN as f64
                - scale_y(
                    element.1 as f64,
                    0.0,
                    height_per_plot as f64,
                    0.0,
                    y_max as f64 + subplot_gap
                )
        );
        bin += x_bin;
    }
    Some(path)
}

/// Add the path elements from Vec<PlotData.path> to their SVG tags.
fn add_all_path_elements(plot_data: Vec<PlotData>, height_subplot: isize, width: i32) -> String {
    let mut all_paths = String::new();
    let length = plot_data.len();
    // gap between subplots; could be added as a parameter.
    let subplot_gap = 25;
    // chromosome labels WRONG as they are being labelled from the top down, rather than the bottom up.
    for (i, row) in plot_data.iter().enumerate() {
        // add text tag here
        all_paths += &format!(
            "<text x='25' y='{}' class='chromosome_label' font-family='monospace'>{}\n↓</text>",
            (i as isize * height_subplot) + subplot_gap - 5 + MARGIN as isize,
            row.id
        );
        // add axis tags here
        all_paths += &format!(
            "<text x='{}' y='{}' class='x_axis_label' font-family='monospace'>{} ↑</text>",
            // 37 aligns the up arrow
            width - 2 * MARGIN - 37,
            (i as isize * height_subplot) + subplot_gap - 5 + MARGIN as isize + height_subplot,
            format_number_to_mb(row.max)
        );
        // reverse the order of the paths!
        all_paths += &format!("<path d='{}' id='{}' class='chromosome_line' stroke='black' fill='none' stroke-width='1' transform='translate(0,{})'/>\n", 
                plot_data[length - i - 1].path,
                plot_data[length - i - 1].id,
                -(i as isize * height_subplot));
    }
    all_paths
}

/// Format [`usize`] to megabase string.
fn format_number_to_mb(n: usize) -> String {
    format!("{:.1}Mb", (n as f64 / 1_000_000_f64))
}

/// The final data structure.
#[derive(Debug, Clone)]
pub struct PlotData {
    /// Chromosome ID.
    pub id: String,
    /// SVG path attribute.
    pub path: String,
    /// Max length of chromosome.
    pub max: usize,
    /// Name of the telomeric repeat (not needed?).
    pub sequence: String,
}

/// Loop through the parsed TSV file and
/// calculate SVG path elements on the fly
/// along with other [`PlotData`] elements.
fn generate_plot_data(
    parsed_tsv: Vec<TelomericRepeatRecord>,
    height: i32,
    width: i32,
    height_per_plot: i32,
) -> Vec<PlotData> {
    // so we can break the loop
    let file_length = parsed_tsv.len();
    // the iteration of the loop
    let mut it = 0usize;
    // a mutable vector to calculate svg path attribute
    let mut path_vec = Vec::new();
    let mut plot_data = Vec::new();
    let mut y_max = 0;

    loop {
        if it == file_length - 1 {
            // there may not be a path element
            // so explicitly make a blank if there is not.
            let path_element = match make_path_element(
                path_vec.clone(),
                path_vec.clone().len(),
                y_max as usize,
                height,
                width,
                height_per_plot,
            ) {
                Some(x) => x,
                None => " ".to_owned(),
            };

            plot_data.push(PlotData {
                id: parsed_tsv[it].id.clone(),
                path: path_element,
                max: parsed_tsv[it].window as usize,
                sequence: parsed_tsv[it].telomeric_repeat.clone(),
            });
            break;
        }

        if parsed_tsv[it].id == parsed_tsv[it + 1].id {
            // calculate y max
            if y_max <= parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number
            {
                y_max = parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number;
            }
            // window (i.e x)
            // forward + reverse counts
            path_vec.push((
                parsed_tsv[it].window,
                parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number,
            ));
            it += 1;
        } else {
            // want to calculate y_max and...
            if y_max <= parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number
            {
                y_max = parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number;
            }
            // the path vector for the last element (seems important for things which occur at the
            // ends of chromosomes right..? DOH)
            path_vec.push((
                parsed_tsv[it].window,
                parsed_tsv[it].forward_repeat_number + parsed_tsv[it].reverse_repeat_number,
            ));
            // calculate the svg path element from path_vec here
            // there may not be a path element
            // so explicitly make a blank if there is not.
            let path_element = match make_path_element(
                path_vec.clone(),
                path_vec.clone().len(),
                y_max as usize,
                height,
                width,
                height_per_plot,
            ) {
                Some(x) => x,
                None => " ".to_owned(),
            };

            plot_data.push(PlotData {
                id: parsed_tsv[it].id.clone(),
                path: path_element,
                max: parsed_tsv[it].window as usize,
                sequence: parsed_tsv[it].telomeric_repeat.clone(),
            });
            path_vec.clear();
            it += 1;
            y_max = 0;
        }
    }
    plot_data
}
