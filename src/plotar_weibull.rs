use crate::weibull::Weibull;
use plotters::prelude::*;
use std::error::Error;

pub fn plotar_weibull(ensaios: &[f64], weibull: &Weibull) -> Result<(), Box<dyn Error>> {
    if ensaios.len() < 2 {
        return Err("São necessários pelo menos dois ensaios.".into());
    }

    let mut dados = ensaios.to_vec();
    dados.sort_by(|a, b| a.total_cmp(b));

    let n = dados.len() as f64;

    let x_min = dados[0] * 0.90;
    let x_max = dados[dados.len() - 1] * 1.10;

    let root = BitMapBackend::new("weibull_ruptura.png", (1100, 750)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Probabilidade de Ruptura - Weibull", ("sans-serif", 30))
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

            Circle::new((sigma, p * 100.0), 5, BLUE.filled())
        }))?
        .label("Dados experimentais")
        .legend(|(x, y)| Circle::new((x, y), 5, BLUE.filled()));

    // Curva teórica de Weibull ajustada.
    let pontos_curva = (0..=300).map(|i| {
        let sigma = x_min + (x_max - x_min) * i as f64 / 300.0;

        (sigma, weibull.cdf(sigma) * 100.0)
    });

    chart
        .draw_series(LineSeries::new(pontos_curva, RED.stroke_width(3)))?
        .label("Weibull ajustada")
        .legend(|(x, y)| PathElement::new(vec![(x - 10, y), (x + 10, y)], RED.stroke_width(3)));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperLeft)
        .draw()?;

    root.present()?;

    println!("\nGráfico salvo em: weibull_ruptura.png");

    Ok(())
}

pub fn plotar_weibull_linearizado(
    ensaios: &[f64],
    weibull: &Weibull,
) -> Result<(), Box<dyn std::error::Error>> {
    if ensaios.len() < 2 {
        return Err("São necessários pelo menos dois ensaios.".into());
    }

    if ensaios.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err("As tensões devem ser finitas e positivas.".into());
    }

    let mut dados = ensaios.to_vec();
    dados.sort_by(|a, b| a.total_cmp(b));

    let n = dados.len() as f64;

    // Transforma a tensão e a probabilidade empírica.
    // Posição de plotagem de Bernard:
    // F_i = (i - 0.3) / (n + 0.4), com i começando em 1.
    let pontos: Vec<(f64, f64)> = dados
        .iter()
        .enumerate()
        .map(|(i, &sigma)| {
            let f = (i as f64 + 1.0 - 0.3) / (n + 0.4);

            let x = sigma.ln();
            let y = (-(1.0 - f).ln()).ln();

            (x, y)
        })
        .collect();

    let x_min = pontos.first().unwrap().0;
    let x_max = pontos.last().unwrap().0;

    let y_min_dados = pontos.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);

    let y_max_dados = pontos.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);

    let margem_y = (y_max_dados - y_min_dados) * 0.15;

    let y_min = y_min_dados - margem_y;
    let y_max = y_max_dados + margem_y;

    let root = BitMapBackend::new("weibull_linearizado.png", (1100, 750)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Gráfico de Probabilidade Weibull Linearizado",
            ("sans-serif", 28),
        )
        .margin(25)
        .x_label_area_size(60)
        .y_label_area_size(80)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc("ln(Tensão de ruptura [MPa])")
        .y_desc("ln(-ln(1 - F))")
        .draw()?;

    // Pontos experimentais transformados.
    chart
        .draw_series(
            pontos
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 5, BLUE.filled())),
        )?
        .label("Dados experimentais")
        .legend(|(x, y)| Circle::new((x, y), 5, BLUE.filled()));

    // Reta teórica calculada com os parâmetros Weibull ajustados:
    // y = k * ln(sigma) - k * ln(lambda).
    let reta = (0..=200).map(|i| {
        let x = x_min + (x_max - x_min) * i as f64 / 200.0;

        let y = weibull.shape * (x - weibull.scale.ln());

        (x, y)
    });

    chart
        .draw_series(LineSeries::new(reta, RED.stroke_width(3)))?
        .label(format!(
            "Weibull ajustada: k = {:.3}, lambda = {:.3} MPa",
            weibull.shape, weibull.scale,
        ))
        .legend(|(x, y)| PathElement::new(vec![(x - 10, y), (x + 10, y)], RED.stroke_width(3)));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperLeft)
        .draw()?;

    root.present()?;

    println!("Gráfico Weibull linearizado salvo em: weibull_linearizado.png");

    Ok(())
}
