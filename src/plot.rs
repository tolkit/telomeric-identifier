pub mod plot {

    // read in csv and make svg plot

    use csv::ReaderBuilder;
    use serde::Deserialize;
    use std::error::Error;
    use std::fs::File;
    use std::io::prelude::*;

    // the CSV records
    #[derive(Debug, Deserialize)]
    pub struct TelomericRepeatRecord {
        pub ID: String,
        pub window: i32,
        pub forward_repeat_number: i32,
        pub reverse_repeat_number: i32,
        pub telomeric_repeat: String,
    }

    fn parse_csv(path: &str) -> Vec<TelomericRepeatRecord> {
        let mut csv_reader = ReaderBuilder::new().from_path(path).unwrap();
        let mut plot_coords_vec = Vec::new();
        for result in csv_reader.deserialize() {
            let record: TelomericRepeatRecord = result.unwrap();
            plot_coords_vec.push(record);
        }
        plot_coords_vec
    }

    fn chromosome_number(parsed_csv: &Vec<TelomericRepeatRecord>) -> usize {
        let mut ids = Vec::new();
        for i in parsed_csv {
            ids.push(&i.ID);
        }
        ids.sort_unstable();
        ids.dedup();
        ids.len()
    }

    fn make_path_element(path_vec: Vec<(i32, i32)>, max: usize, height: i32, width: i32) -> String {
        // somehow going to have to scale the lines to fit on all graphs
        const margin: i32 = 40;
        let mut path = String::new();
        let width_incl_margins = width - margin;
        let heigh_incl_margins = height - margin;

        let x_bin: f64 = width_incl_margins as f64 / max as f64;

        // add the first move to point
        path += &format!(
            "M{},{}",
            (margin / 2) as f64 + x_bin,
            height - path_vec[0].1
        );

        let mut bin: f64 = 0.0;
        for element in path_vec.iter().skip(1) {
            path += &format!("L{},{}", bin + (margin as f64 / 2.0), height - element.1);
            bin += x_bin;
        }
        path
    }

    fn add_all_path_elements(plot_data: Vec<PlotData>) -> String {
        let mut all_paths = String::new();
        for (i, row) in plot_data.iter().enumerate() {
            all_paths += &format!("<path d='{}' stroke='black' fill='none' stroke-width='1' transform='translate(0,{})'/>\n", row.path, i as isize * -150);
        }
        all_paths
    }

    #[derive(Debug, Clone)]
    pub struct PlotData {
        // chromosome
        pub id: String,
        // svg path attribute
        pub path: String,
        // max length of chromosome
        pub max: usize,
        // name of the telomeric repeat (not needed?)
        pub sequence: String,
    }

    fn generate_plot_data(
        parsed_csv: Vec<TelomericRepeatRecord>,
        height: i32,
        width: i32,
    ) -> Vec<PlotData> {
        //-> Vec<PlotData>
        // so we can break the loop
        let file_length = parsed_csv.len();
        // the iteration of the loop
        let mut it = 0usize;
        // a mutable vector to calculate svg path attribute
        let mut path_vec = Vec::new();
        let mut plot_data = Vec::new();

        loop {
            if it == file_length - 1 {
                plot_data.push(PlotData {
                    id: parsed_csv[it].ID.clone(),
                    path: make_path_element(
                        path_vec.clone(),
                        path_vec.clone().len(),
                        height,
                        width,
                    ),
                    max: parsed_csv[it].window as usize,
                    sequence: parsed_csv[it].telomeric_repeat.clone(),
                });
                break;
            }

            if parsed_csv[it].ID == parsed_csv[it + 1].ID {
                // window (i.e x)
                // forward + reverse counts
                path_vec.push((
                    parsed_csv[it].window,
                    parsed_csv[it].forward_repeat_number + parsed_csv[it].reverse_repeat_number,
                ));
                it += 1;
            } else {
                // calculate the svg path element from path_vec here

                plot_data.push(PlotData {
                    id: parsed_csv[it].ID.clone(),
                    path: make_path_element(
                        path_vec.clone(),
                        path_vec.clone().len(),
                        height,
                        width,
                    ),
                    max: parsed_csv[it].window as usize,
                    sequence: parsed_csv[it].telomeric_repeat.clone(),
                });
                path_vec.clear();
                it += 1;
            }
        }

        plot_data
    }

    pub fn plot(matches: &clap::ArgMatches) -> std::io::Result<()> {
        let csv = matches.value_of("csv").unwrap();

        let parsed_csv = parse_csv(csv);

        // calculate this instead from the plot data, after
        // filtering on chromosomes < XX Mb
        let number_of_chromosomes = chromosome_number(&parsed_csv);

        let width = 1000;
        // allow height 150 per chromosome?
        let height: i32 = 150 * number_of_chromosomes as i32;

        let plot_data = generate_plot_data(parsed_csv, height, width);

        for i in plot_data.clone() {
            println!("{:?}", i);
        }

        let out_filename = format!("test.svg");

        let mut svg_file = File::create(out_filename)?;

        let svg = format!(
            "<?xml version='1.0' encoding='UTF-8'  standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'> \
                 {} \
                 </svg>",
            width,
            height,
            add_all_path_elements(plot_data)
        );

        svg_file.write_all(svg.as_bytes()).expect("unable to write");

        Ok(())
    }
}
