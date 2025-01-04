use ark_ff::{AdditiveGroup, Field};
use crate::{constants::CircomPoseidonConstants, field::Fr};

pub struct CircomPoseidon {
    constants: CircomPoseidonConstants,
    state: [Fr; 3],
}

impl CircomPoseidon {
    pub fn hash(x: Fr, y: Fr) -> Fr {
        let constants = CircomPoseidonConstants::default();
        let initial_state = [Fr::default(), x, y];
        let initial = CircomPoseidon {
            constants,
            state: initial_state,
        };

        let state_after_first_ark = initial.apply_ark(0);
        let processed_state = (0..(8 / 2 - 1)).fold(state_after_first_ark, |current, r| {
            let next_state = current.apply_sigma();
            let after_ark = next_state.apply_ark((r + 1) * 3);
            after_ark.apply_mix_m()
        });

        let after_partial_round = processed_state
            .apply_sigma()
            .apply_ark((8 / 2) * 3)
            .apply_mix_p();

        let after_full_rounds = (0..57).fold(after_partial_round, |current, r| {
            let updated_state = [
                CircomPoseidon::sigma(current.state[0]) + current.constants.c[(8 / 2 + 1) * 3 + r],
                current.state[1],
                current.state[2],
            ];
            let mixed_state = CircomPoseidon {
                state: updated_state,
                constants: current.constants,
            };
            mixed_state.apply_mix_s(r)
        });

        let final_processed_state = (0..(8 / 2 - 1)).fold(after_full_rounds, |current, r| {
            let next_state = current.apply_sigma();
            let after_ark = next_state.apply_ark((8 / 2 + 1) * 3 + 57 + r * 3);
            after_ark.apply_mix_m()
        });

        final_processed_state.apply_sigma().mix_last(0)
    }

    fn sigma(x: Fr) -> Fr {
        x.pow([5])
    }

    fn apply_ark(self, r: usize) -> Self {
        let state = self
            .state
            .iter()
            .enumerate()
            .map(|(i, &x)| x + self.constants.c[i + r])
            .collect::<Vec<Fr>>()
            .try_into()
            .unwrap();

        Self { state, ..self }
    }

    fn apply_sigma(self) -> Self {
        let state = self
            .state
            .iter()
            .map(|&x| CircomPoseidon::sigma(x))
            .collect::<Vec<Fr>>()
            .try_into()
            .unwrap();

        Self { state, ..self }
    }

    fn apply_mix_m(self) -> Self {
        let state = (0..self.state.len())
            .map(|j| {
                self.state
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| x * self.constants.m[i][j])
                    .fold(Fr::ZERO, |acc, val| acc + val)
            })
            .collect::<Vec<Fr>>()
            .try_into()
            .unwrap();

        Self { state, ..self }
    }

    fn apply_mix_p(self) -> Self {
        let state = (0..self.state.len())
            .map(|j| {
                self.state
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| x * self.constants.p[i][j])
                    .fold(Fr::ZERO, |acc, val| acc + val)
            })
            .collect::<Vec<Fr>>()
            .try_into()
            .unwrap();

        Self { state, ..self }
    }

    fn apply_mix_s(self, r: usize) -> Self {
        let state = [
            self.state
                .iter()
                .cloned()
                .enumerate()
                .map(|(i, x)| x * self.constants.s[(3 * 2 - 1) * r + i])
                .fold(Fr::ZERO, |acc, val| acc + val),
        ]
        .into_iter()
        .chain((1..self.state.len()).map(|i| {
            self.state[i] + self.state[0] * self.constants.s[(3 * 2 - 1) * r + 3 + i - 1]
        }))
        .collect::<Vec<Fr>>()
        .try_into()
        .unwrap();

        Self { state, ..self }
    }

    fn mix_last(self, s: usize) -> Fr {
        self.state
            .iter()
            .enumerate()
            .map(|(i, &x)| x * self.constants.m[i][s])
            .fold(Fr::ZERO, |acc, val| acc + val)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use ark_ff::{BigInt, PrimeField};

    use super::*;

    #[test]
    fn test_sigma() {
        let x = Fr::from_bigint(
            BigInt::from_str(
                "10498446744258291650735421565632985639701981041287181768865869454200878131887",
            )
            .unwrap(),
        )
        .unwrap();
        let result = CircomPoseidon::sigma(x);
        assert_eq!(
            result.into_bigint(),
            BigInt::from_str(
                "12482121694093354782765194089691070098368127388830956394165914815968551446065",
            )
            .unwrap()
        );
    }

    #[test]
    fn test_mix_m() {
        let constants = CircomPoseidonConstants::default();
        let s = CircomPoseidon {
            constants,
            state: [
                Fr::from_bigint(
                    BigInt::from_str(
                        "18748214218909041356337553469301793477858404099161101255670448756216508301025",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "14001090959915885276730193140615367824457344208709298148449921602972961864271",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "21039621298411050525972067638905940728497671706374440393718856029109245780568",
                    )
                    .unwrap(),
                )
                .unwrap(),
            ]
        };

        let result = s.apply_mix_m();

        let expected_result = [
            Fr::from_bigint(
                BigInt::from_str(
                    "19212249459613001635322677932044571594832098024684187537767379651838537722242",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "13346432103628585998281534921408530751039791002079061852835106438825147371757",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "470145441107885699228147589021085756242352387945133942451294769220019691807",
                )
                .unwrap(),
            )
            .unwrap(),
        ];
        assert_eq!(result.state, expected_result);
    }

    #[test]
    fn test_mix_p() {
        let constants = CircomPoseidonConstants::default();
        let s = CircomPoseidon {
            constants,
            state: [
                Fr::from_bigint(
                    BigInt::from_str(
                        "11054539610865716959112730216921823956884402188047559139840762443903078964807",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "889496775026312796519119247885544826173520850299645133418342221246834572377",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "20327605851775715928805880379802842509685187839623124393912839998555298548830",
                    )
                    .unwrap(),
                )
                .unwrap(),
            ]
        };

        let result = s.apply_mix_p();

        let expected_result = [
            Fr::from_bigint(
                BigInt::from_str(
                    "10957445925160270044693914644787882019446927193081715019239813790612436710990",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "15617724767560381146854323530626642274265602408002782166490987198203153486720",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "8188576910493806248734793780989996031043527802574000617472436041968065581651",
                )
                .unwrap(),
            )
            .unwrap(),
        ];
        assert_eq!(result.state, expected_result);
    }

    #[test]
    fn test_mix_s() {
        let constants = CircomPoseidonConstants::default();
        let s = CircomPoseidon {
            constants,
            state: [
                Fr::from_bigint(
                    BigInt::from_str(
                        "1073104925714253051670422979696287993089561378609784836524058216953452518800",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "18698081889797527287892615152858593919908252735314388801102459459095355586529",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "964734953777785051088731259326514402767379753056823320461487398194316968472",
                    )
                    .unwrap(),
                )
                .unwrap(),
            ]
        };

        let result = s.apply_mix_s(0);

        let expected_result = [
            Fr::from_bigint(
                BigInt::from_str(
                    "4413576681639644708927435504116008819690219979742107806543982128344816486059",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "18047144526083855362077562050797531820746747110365503012268045400381426067118",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "14512113857556374279800803178193270419547881978169561761995632044104001997174",
                )
                .unwrap(),
            )
            .unwrap(),
        ];
        assert_eq!(result.state, expected_result);
    }

    #[test]
    fn test_mix_last() {
        let constants = CircomPoseidonConstants::default();
        let s = CircomPoseidon {
            constants,
            state: [
                Fr::from_bigint(
                    BigInt::from_str(
                        "14903581272446463225204524721363305270359053220610820269937666541151177015349",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "17993715061665243415477544347824042376858770491012817421933497650720324679864",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "15253365088956611786925658848458340812843569370338269931949433692506640669343",
                    )
                    .unwrap(),
                )
                .unwrap(),
            ]
        };

        let result = s.mix_last(0);

        let expected_result = Fr::from_bigint(
            BigInt::from_str(
                "10776602241119476978803538453188198564682373569005749437912787513019500133942",
            )
            .unwrap(),
        )
        .unwrap();
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_ark() {
        let constants = CircomPoseidonConstants::default();
        let s = CircomPoseidon {
            constants,
            state: [
                Fr::from_bigint(
                    BigInt::from_str(
                        "13113380046921429357968206206413495707844213025058902153087258344566081685705",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "15502287019891195293542473848301291711293573680196308374326156854444666425933",
                    )
                    .unwrap(),
                )
                .unwrap(),
                Fr::from_bigint(
                    BigInt::from_str(
                        "960782068423008155923835610511860070078873921748828059626034504954890453563",
                    )
                    .unwrap(),
                )
                .unwrap(),
            ]
        };

        let result = s.apply_ark(0);

        let expected_result = [
            Fr::from_bigint(
                BigInt::from_str(
                    "19858578037131633956342249035175485304147089324604866555944670074438212720439",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "15928568697651131885563790657366470529141658358875818949042050993134916565681",
                )
                .unwrap(),
            )
            .unwrap(),
            Fr::from_bigint(
                BigInt::from_str(
                    "4974970831339591754812778277936825500366371746378485279433975965182263031344",
                )
                .unwrap(),
            )
            .unwrap(),
        ];
        assert_eq!(result.state, expected_result);
    }

    #[test]
    fn test_poseidon() {
        let x = Fr::from_bigint(
            BigInt::from_str(
                "4919343374109933223627495853328117416310439472367144134356200729566675210393",
            )
            .unwrap(),
        )
        .unwrap();

        let y = Fr::from_bigint(
            BigInt::from_str(
                "5253109015430049131200882288700211126182357120559979176158196848667194692355",
            )
            .unwrap(),
        )
        .unwrap();

        let result = CircomPoseidon::hash(x, y);

        let expected_result = Fr::from_bigint(
            BigInt::from_str(
                "4804190683062966564094284552215482492423412176373545691806986132546823080430",
            )
            .unwrap(),
        )
        .unwrap();
        assert_eq!(result, expected_result);
    }
}
