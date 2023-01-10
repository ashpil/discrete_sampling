use exr::prelude::*;

use pdf_maker::Data2D;
use pdf_maker::ContinuousAlias2D;
use pdf_maker::Inversion2D;
use pdf_maker::Hierarchical2D;
use pdf_maker::distribution::Continuous2D;

fn luminance([r, g, b]: [f32; 3]) -> f32 {
    r * 0.2126 + g * 0.7152 + b * 0.0722
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: warping <filename.exr>");
        std::process::exit(1);
    }

    let density_image = read_first_rgba_layer_from_file(args[1].clone(), |resolution, _| {
        let width = resolution.width();
        let height = resolution.height();
        Data2D::new_same(width, height, 0.0)
    }, |buffer, pos, (r, g, b, _): (f32, f32, f32, f32)| {
        buffer[pos.y()][pos.x()] = luminance([r, g, b]);
    }).unwrap().layer_data.channel_data.pixels;


    fn visualize_warping<D: Continuous2D<Weight=f32>>(out_name: &str, weights: &Data2D<f32>) {
        let sampler = D::build(weights);
        let warping = sampler.visualize_warping(16);
        write_rgb_file(out_name, warping.width(), warping.height(), |x, y| {
            let p = warping[y][x];
            (p[0], p[1], p[2])
        }).unwrap();
        println!("Wrote {}", out_name);
    }

    visualize_warping::<Inversion2D<f32>>("inversion_warping.exr", &density_image);
    visualize_warping::<ContinuousAlias2D<f32>>("alias_warping.exr", &density_image);
    visualize_warping::<Hierarchical2D<f32>>("hierarchical_warping.exr", &density_image);
}

