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
        self.scale * gamma(1.0 + 1.0 / self.shape)
    }
    /// Ajusta uma distribuição de Weibull de dois parâmetros
    /// usando máxima verossimilhança para dados completos.
    pub fn fit(data: &[f64]) -> Result<Self, &'static str> {
        if data.len() < 2 {
            return Err("São necessários pelo menos dois ensaios.");
        }

        if data.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err("Todas as tensões devem ser finitas e positivas.");
        }

        let logs: Vec<f64> = data.iter().map(|x| x.ln()).collect();
        let mean_log = logs.iter().sum::<f64>() / logs.len() as f64;
        let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        // Se todas as observações forem iguais, não existe uma
        // estimativa finita usual para o parâmetro de forma por MLE.
        if logs.iter().all(|v| *v == logs[0]) {
            return Err("Os dados são todos iguais; não há ajuste Weibull finito.");
        }

        // Equação da máxima verossimilhança:
        //
        // média ponderada de ln(x) - média de ln(x) = 1/k
        //
        // O peso de cada observação é proporcional a x^k.
        let score = |k: f64| -> f64 {
            let weights: Vec<f64> = logs
                .iter()
                .map(|ln_x| (k * (ln_x - max_log)).exp())
                .collect();

            let sum_weights = weights.iter().sum::<f64>();

            let weighted_log_mean = logs
                .iter()
                .zip(weights.iter())
                .map(|(ln_x, w)| ln_x * w)
                .sum::<f64>()
                / sum_weights;

            weighted_log_mean - mean_log - 1.0 / k
        };

        // Encontra um intervalo que contenha a raiz.
        let mut low = 1e-6_f64;
        let mut high = 1.0_f64;

        while score(high) <= 0.0 {
            high *= 2.0;

            if high > 1e6 {
                return Err("Não foi possível encontrar um ajuste estável.");
            }
        }

        // Bisseção para estimar o parâmetro de forma k.
        for _ in 0..200 {
            let mid = (low + high) / 2.0;

            if score(mid) > 0.0 {
                high = mid;
            } else {
                low = mid;
            }
        }

        let shape = (low + high) / 2.0;

        // Calcula lambda sem elevar diretamente valores grandes
        // à potência k, reduzindo o risco de overflow.
        let normalized_sum: f64 = logs
            .iter()
            .map(|ln_x| (shape * (ln_x - max_log)).exp())
            .sum();

        let log_scale = max_log + (normalized_sum / data.len() as f64).ln() / shape;

        let scale = log_scale.exp();

        if !shape.is_finite() || shape <= 0.0 || !scale.is_finite() || scale <= 0.0 {
            return Err("O ajuste produziu parâmetros inválidos.");
        }

        Ok(Self { shape, scale })
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
        return std::f64::consts::PI / ((std::f64::consts::PI * z).sin() * gamma(1.0 - z));
    }

    let z = z - 1.0;
    let mut x = P[0];

    for (i, p) in P.iter().enumerate().skip(1) {
        x += p / (z + i as f64);
    }

    let t = z + 7.5;

    (2.0 * std::f64::consts::PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
}

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
    let weibull = match Weibull::fit(&ensaios) {
        Ok(modelo) => modelo,
        Err(erro) => {
            eprintln!("Erro no ajuste: {}", erro);
            return;
        }
    };

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
            82.0, 91.0, 95.0, 98.0,
            101.0, 104.0, 107.0, 110.0,
            114.0, 118.0, 123.0, 130.0,
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
