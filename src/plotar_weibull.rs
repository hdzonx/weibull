use plotters::prelude::*;
use std::error::Error;
use crate::weibull::Weibull;

pub fn plotar_weibull(
    ensaios: &[f64],
    weibull: &Weibull,
) -> Result<(), Box<dyn Error>> {
    if ensaios.len() < 2 {
        return Err("São necessários pelo menos dois ensaios.".into());
    }

    let mut dados = ensaios.to_vec();
    dados.sort_by(|a, b| a.total_cmp(b));

    let n = dados.len() as f64;

    let x_min = dados[0] * 0.90;
    let x_max = dados[dados.len() - 1] * 1.10;

    let root = BitMapBackend::new(
        "weibull_ruptura.png",
        (1100, 750),
    )
    .into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Probabilidade de Ruptura - Weibull",
            ("sans-serif", 30),
        )
        .margin(25)
        .x_label_area_size(55)
        .y_label_area_size(75)
        .build_cartesian_2d(x_min..x_max, 0.0..100.0)?;

    chart
        .configure_mesh()
        .x_desc("Tensao de ruptura (MPa)")
        .y_desc("Probabilidade acumulada de ruptura (%)")
        .y_label_formatter(&|v| format!("{:.0}%", v))
        .draw()?;

    // Pontos experimentais usando a posição de plotagem
    // de Bernard: (i - 0.3) / (n + 0.4).
    chart
        .draw_series(dados.iter().enumerate().map(|(i, &sigma)| {
            let p = (i as f64 + 1.0 - 0.3) / (n + 0.4);

            Circle::new(
                (sigma, p * 100.0),
                5,
                BLUE.filled(),
            )
        }))?
        .label("Dados experimentais")
        .legend(|(x, y)| {
            Circle::new((x, y), 5, BLUE.filled())
        });

    // Curva teórica de Weibull ajustada.
    let pontos_curva = (0..=300).map(|i| {
        let sigma = x_min
            + (x_max - x_min) * i as f64 / 300.0;

        (sigma, weibull.cdf(sigma) * 100.0)
    });

    chart
        .draw_series(LineSeries::new(
            pontos_curva,
            RED.stroke_width(3),
        ))?
        .label("Weibull ajustada")
        .legend(|(x, y)| {
            PathElement::new(
                vec![(x - 10, y), (x + 10, y)],
                RED.stroke_width(3),
            )
        });

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperLeft)
        .draw()?;

    root.present()?;

    println!(
        "\nGráfico salvo em: weibull_ruptura.png"
    );

    Ok(())
}