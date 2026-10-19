// Distribuição de Weibull de dois parâmetros.
#[derive(Debug, Clone, Copy)]
pub struct Weibull {
    pub shape: f64, // k: parâmetro de forma
    pub scale: f64, // lambda: parâmetro de escala
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
