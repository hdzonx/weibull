use std::f64::consts::E;
pub mod plotar_weibull;
pub mod weibull;

fn main() {
    // Tensões de ruptura medidas nos corpos de prova (MPa).
    let ensaios: Vec<f64> = vec![
        82.0, 91.0, 95.0, 98.0, 101.0, 104.0, 107.0, 110.0, 114.0, 118.0, 123.0, 130.0,
    ];

    println!("ANÁLISE DE RESISTÊNCIA DOS CORPOS DE PROVA");
    println!("==========================================");

    println!("Quantidade de ensaios: {}", ensaios.len());

    let media = ensaios.iter().sum::<f64>() / ensaios.len() as f64;

    println!("Tensão média observada: {:.2} MPa", media);

    // Ajuste automático da distribuição de Weibull.
    let weibull = match weibull::Weibull::fit(&ensaios) {
        Ok(modelo) => modelo,
        Err(erro) => {
            eprintln!("Erro no ajuste: {}", erro);
            return;
        }
    };

    if let Err(erro) = plotar_weibull::plotar_weibull(&ensaios, &weibull) {
        eprintln!("Erro ao gerar gráfico: {}", erro);
    }

    if let Err(erro) = plotar_weibull::plotar_weibull_linearizado(&ensaios, &weibull) {
        eprintln!("Erro no gráfico linearizado: {}", erro);
    }

    println!("\nPARÂMETROS ESTIMADOS POR MÁXIMA VEROSSIMILHANÇA");
    println!("k (forma): {:.6}", weibull.shape);
    println!("lambda (escala): {:.6} MPa", weibull.scale);

    println!("\nPROBABILIDADES ESTIMADAS PELO MODELO");
    println!(
        "{:>14} {:>18} {:>18}",
        "Tensão (MPa)", "Ruptura acumulada", "Sobrevivência"
    );

    for sigma in [80.0, 90.0, 100.0, 110.0, 120.0, 130.0] {
        println!(
            "{:>14.1} {:>17.3}% {:>17.3}%",
            sigma,
            weibull.cdf(sigma) * 100.0,
            weibull.survival(sigma) * 100.0
        );
    }

    // Resistência associada a 5% de ruptura acumulada.
    let p: f64 = 0.05;
    let sigma_5 = weibull.scale * (-(1.0_f64 - p).ln()).powf(1.0 / weibull.shape);

    println!("\nTensão associada a 5% de ruptura: {:.2} MPa", sigma_5);

    // Resistência associada a 10% de ruptura acumulada.
    let p10: f64 = 0.10;
    let sigma_10 = weibull.scale * (-(1.0_f64 - p10).ln()).powf(1.0 / weibull.shape);

    println!("Tensão associada a 10% de ruptura: {:.2} MPa", sigma_10);
}

#[cfg(test)]
mod tests {
    use super::*;
    use weibull::Weibull;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "Esperado {b}, obtido {a}");
    }

    #[test]
    fn test_parametros_validos() {
        assert!(Weibull::new(2.0, 1000.0).is_ok());
    }

    #[test]
    fn test_parametros_invalidos() {
        assert!(Weibull::new(0.0, 1000.0).is_err());
        assert!(Weibull::new(-1.0, 1000.0).is_err());
        assert!(Weibull::new(f64::NAN, 1000.0).is_err());
        assert!(Weibull::new(2.0, 0.0).is_err());
        assert!(Weibull::new(2.0, f64::INFINITY).is_err());
    }

    #[test]
    fn test_cdf_em_zero_e_negativo() {
        let w = Weibull::new(2.0, 1000.0).unwrap();

        close(w.cdf(0.0), 0.0, 1e-12);
        close(w.cdf(-10.0), 0.0, 1e-12);
    }

    #[test]
    fn test_cdf_e_sobrevivencia() {
        let w = Weibull::new(2.0, 1000.0).unwrap();
        let x = 500.0;

        let esperado = 1.0 - (-0.25_f64).exp();

        close(w.cdf(x), esperado, 1e-12);
        close(w.survival(x), (-0.25_f64).exp(), 1e-12);
        close(w.cdf(x) + w.survival(x), 1.0, 1e-12);
    }

    #[test]
    fn test_pdf_exponencial() {
        // Weibull com k = 1 equivale à distribuição exponencial.
        let w = Weibull::new(1.0, 2.0).unwrap();

        close(w.pdf(0.0), 0.5, 1e-12);
        close(w.pdf(2.0), 0.5 / E, 1e-12);
        close(w.cdf(2.0), 1.0 - 1.0 / E, 1e-12);
    }

    #[test]
    fn test_pdf_na_origem() {
        let k_menor = Weibull::new(0.5, 2.0).unwrap();
        let k_igual = Weibull::new(1.0, 2.0).unwrap();
        let k_maior = Weibull::new(2.0, 2.0).unwrap();

        assert_eq!(k_menor.pdf(0.0), f64::INFINITY);
        close(k_igual.pdf(0.0), 0.5, 1e-12);
        close(k_maior.pdf(0.0), 0.0, 1e-12);
    }

    #[test]
    fn test_media_teorica() {
        // Para k = 1, a média é igual a lambda.
        let w = Weibull::new(1.0, 10.0).unwrap();

        close(w.mean(), 10.0, 1e-10);
    }

    #[test]
    fn test_limites_de_probabilidade() {
        let w = Weibull::new(2.0, 1000.0).unwrap();

        for x in [0.0, 100.0, 500.0, 1000.0, 5000.0] {
            let cdf = w.cdf(x);
            assert!((0.0..=1.0).contains(&cdf));
        }

        assert_eq!(w.pdf(-1.0), 0.0);
        assert_eq!(w.survival(-1.0), 1.0);
    }

    #[test]
    fn ajusta_dados_positivos() {
        let dados = [
            82.0, 91.0, 95.0, 98.0, 101.0, 104.0, 107.0, 110.0, 114.0, 118.0, 123.0, 130.0,
        ];

        let w = Weibull::fit(&dados).unwrap();

        assert!(w.shape.is_finite() && w.shape > 0.0);
        assert!(w.scale.is_finite() && w.scale > 0.0);

        // A CDF deve ser crescente e as probabilidades
        // devem permanecer entre zero e um.
        let mut anterior = 0.0;

        for x in [80.0, 90.0, 100.0, 110.0, 120.0, 130.0] {
            let p = w.cdf(x);

            assert!((0.0..=1.0).contains(&p));
            assert!(p >= anterior);

            anterior = p;
        }
    }

    #[test]
    fn rejeita_dados_invalidos() {
        assert!(Weibull::fit(&[]).is_err());
        assert!(Weibull::fit(&[100.0]).is_err());
        assert!(Weibull::fit(&[80.0, 0.0, 100.0]).is_err());
        assert!(Weibull::fit(&[80.0, -10.0]).is_err());
        assert!(Weibull::fit(&[80.0, f64::NAN]).is_err());
        assert!(Weibull::fit(&[80.0, f64::INFINITY]).is_err());
    }

    #[test]
    fn rejeita_observacoes_identicas() {
        assert!(Weibull::fit(&[100.0, 100.0, 100.0]).is_err());
    }

    #[test]
    fn escala_fica_proxima_da_media_quando_k_e_um() {
        // Para dados exponenciais, a Weibull com k = 1
        // equivale a uma distribuição exponencial.
        // Este teste verifica a consistência das funções
        // de distribuição para um modelo de k = 1.
        let w = Weibull::new(1.0, 100.0).unwrap();

        assert!((w.cdf(100.0) - (1.0 - (-1.0_f64).exp())).abs() < 1e-12);
    }
}
