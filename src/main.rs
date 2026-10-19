
use std::f64::consts::E;

/// Distribuição de Weibull de dois parâmetros.
#[derive(Debug, Clone, Copy)]
pub struct Weibull {
    shape: f64, // k: parâmetro de forma
    scale: f64, // lambda: parâmetro de escala
}

impl Weibull {
    /// Cria uma distribuição de Weibull.
    pub fn new(shape: f64, scale: f64) -> Result<Self, &'static str> {
        if !shape.is_finite() || shape <= 0.0 {
            return Err("O parâmetro shape deve ser finito e positivo.");
        }

        if !scale.is_finite() || scale <= 0.0 {
            return Err("O parâmetro scale deve ser finito e positivo.");
        }

        Ok(Self { shape, scale })
    }

    /// Função densidade de probabilidade (PDF).
    pub fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }

        let k = self.shape;
        let lambda = self.scale;

        if x == 0.0 {
            return if k < 1.0 {
                f64::INFINITY
            } else if k == 1.0 {
                1.0 / lambda
            } else {
                0.0
            };
        }

        let z = x / lambda;
        (k / lambda) * z.powf(k - 1.0) * (-(z.powf(k))).exp()
    }

    /// Função de distribuição acumulada (CDF).
    pub fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }

        let z = x / self.scale;
        -(-(z.powf(self.shape))).exp_m1()
    }

    /// Função de sobrevivência: P(T > x).
    pub fn survival(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 1.0;
        }

        let z = x / self.scale;
        (-(z.powf(self.shape))).exp()
    }

    /// Taxa de falha instantânea (função de risco).
    pub fn hazard(&self, x: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }

        self.pdf(x) / self.survival(x)
    }

    /// Média teórica da distribuição.
    pub fn mean(&self) -> f64 {
        self.scale
            * gamma(1.0 + 1.0 / self.shape)
    }
}

/// Aproximação da função Gamma pela fórmula de Lanczos.
fn gamma(z: f64) -> f64 {
    const P: [f64; 9] = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];

    if z < 0.5 {
        return std::f64::consts::PI
            / ((std::f64::consts::PI * z).sin()
                * gamma(1.0 - z));
    }

    let z = z - 1.0;
    let mut x = P[0];

    for (i, p) in P.iter().enumerate().skip(1) {
        x += p / (z + i as f64);
    }

    let t = z + 7.5;

    (2.0 * std::f64::consts::PI).sqrt()
        * t.powf(z + 0.5)
        * (-t).exp()
        * x
}


fn main() {
    // Resultados experimentais hipotéticos:
    // tensões de ruptura em MPa.
    let ensaios: Vec<f64> = vec![
        82.0, 91.0, 95.0, 98.0,
        101.0, 104.0, 107.0, 110.0,
        114.0, 118.0, 123.0, 130.0,
    ];

    println!("ANÁLISE DE RUPTURA DOS CORPOS DE PROVA");
    println!("======================================");

    println!("Quantidade de ensaios: {}", ensaios.len());

    let media = ensaios.iter().sum::<f64>()
        / ensaios.len() as f64;

    println!("Tensão média de ruptura: {:.2} MPa", media);

    let menor = ensaios.iter()
        .copied()
        .fold(f64::INFINITY, f64::min);

    let maior = ensaios.iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    println!("Menor tensão de ruptura: {:.2} MPa", menor);
    println!("Maior tensão de ruptura: {:.2} MPa", maior);

    // Parâmetros ilustrativos, não ajustados aos ensaios.
    let shape = 8.0;
    let scale = 110.0;

    let weibull = Weibull::new(shape, scale).unwrap();

    println!("\nMODELO DE WEIBULL ILUSTRATIVO");
    println!("Forma k: {:.2}", shape);
    println!("Escala lambda: {:.2} MPa", scale);

    println!(
        "\n{:>12} {:>18} {:>20}",
        "Tensão (MPa)",
        "Ruptura acumulada",
        "Sobrevivência"
    );

    for sigma in [80.0, 90.0, 100.0, 110.0, 120.0, 130.0] {
        println!(
            "{:>12.1} {:>17.2}% {:>19.2}%",
            sigma,
            weibull.cdf(sigma) * 100.0,
            weibull.survival(sigma) * 100.0
        );
    }

    // Tensão correspondente a uma probabilidade
    // acumulada de ruptura de 5%.
    let p:f64 = 0.05;
    let sigma_5 = scale * (-(1.0 - p).ln()).powf(1.0 / shape);

    println!(
        "\nTensão correspondente a 5% de ruptura: {:.2} MPa",
        sigma_5
    );
}


#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!(
            (a - b).abs() < tol,
            "Esperado {b}, obtido {a}"
        );
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
}
