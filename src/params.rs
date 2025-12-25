use fp2::traits::Fp as _;
use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX};

use crate::bn::BigNum;

pub mod poke_i {
    use super::*;
    use crate::{
        fields::{PokeFieldI, PokeFieldIBase},
        poke::PublicParams,
    };

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const P_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        156, 1, 55, 113, 211, 191, 79, 224, 97, 54, 107, 37, 254, 167, 210, 138, 199, 125, 108,
        159, 62, 27, 61, 12, 176, 93, 127, 206, 236, 40, 77, 235, 18, 81, 163, 191, 61, 216, 30,
        105, 141, 244, 112, 38, 122, 199, 207, 251, 158, 170, 70, 187, 238, 73,
    ];
    const Q_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const Q_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        99, 254, 200, 142, 44, 64, 176, 31, 158, 201, 148, 218, 1, 88, 45, 117, 138, 22, 138, 198,
        255, 79, 247, 48, 32, 60, 177, 74, 18, 196, 166, 147, 100, 195, 53, 43, 189, 187, 224, 237,
        23, 35, 156, 158, 249, 126, 66, 202, 241, 91, 209, 56, 29, 32,
    ];
    const PQ_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        84, 5, 39, 47, 28, 26, 72, 53, 222, 25, 244, 169, 18, 138, 123, 250, 227, 91, 135, 191,
        182, 168, 208, 156, 231, 66, 10, 171, 57, 90, 207, 9, 222, 195, 240, 102, 7, 222, 148, 122,
        208, 175, 249, 130, 55, 245, 12, 92, 175, 174, 252, 231, 208, 15,
    ];
    const PQ_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [0; PokeFieldIBase::ENCODED_LENGTH];
    const P_X: PokeFieldI = PokeFieldI::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldI = PokeFieldI::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const R_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        242, 78, 54, 243, 244, 158, 40, 209, 213, 36, 144, 132, 126, 115, 146, 252, 115, 95, 79,
        49, 121, 90, 228, 120, 114, 82, 233, 129, 214, 22, 113, 22, 116, 81, 115, 222, 238, 180,
        157, 29, 159, 205, 134, 216, 253, 65, 214, 79, 148, 149, 147, 24, 195, 7,
    ];
    const S_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const S_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        13, 177, 201, 12, 11, 97, 215, 46, 42, 219, 111, 123, 129, 140, 109, 3, 222, 52, 167, 52,
        197, 16, 80, 196, 93, 71, 71, 151, 40, 214, 130, 104, 3, 195, 101, 12, 12, 223, 97, 57, 6,
        74, 134, 236, 117, 4, 60, 118, 252, 112, 132, 219, 72, 98,
    ];
    const RS_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        214, 198, 30, 107, 72, 179, 73, 69, 35, 50, 116, 38, 27, 143, 85, 161, 55, 176, 109, 176,
        64, 247, 227, 127, 52, 115, 253, 72, 217, 177, 78, 213, 224, 192, 75, 192, 253, 45, 130,
        177, 170, 220, 184, 89, 185, 137, 120, 89, 231, 163, 80, 255, 92, 95,
    ];
    const RS_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [0; PokeFieldIBase::ENCODED_LENGTH];
    const R_X: PokeFieldI = PokeFieldI::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldI = PokeFieldI::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldI = PokeFieldI::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Construct basis points for the 5^c-torsion on E_0
    const X_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const X_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        211, 24, 39, 150, 242, 9, 42, 168, 202, 224, 82, 61, 102, 182, 231, 124, 213, 107, 144, 72,
        84, 15, 181, 210, 65, 156, 234, 60, 141, 56, 253, 222, 254, 41, 3, 136, 237, 101, 182, 89,
        189, 117, 17, 158, 8, 209, 192, 197, 185, 255, 80, 133, 107, 67,
    ];
    const Y_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const Y_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        44, 231, 216, 105, 13, 246, 213, 87, 53, 31, 173, 194, 153, 73, 24, 131, 124, 40, 102, 29,
        234, 91, 127, 106, 142, 253, 69, 220, 113, 180, 246, 159, 120, 234, 213, 98, 13, 46, 73,
        253, 231, 161, 251, 38, 107, 117, 81, 0, 215, 6, 199, 110, 160, 38,
    ];
    const XY_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        7, 127, 166, 170, 67, 68, 236, 216, 13, 61, 70, 222, 190, 115, 147, 244, 207, 140, 116,
        141, 195, 61, 63, 202, 239, 236, 93, 15, 92, 242, 111, 151, 53, 67, 144, 196, 218, 77, 91,
        160, 138, 225, 199, 32, 138, 40, 48, 229, 231, 107, 234, 178, 109, 61,
    ];
    const XY_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [0; PokeFieldIBase::ENCODED_LENGTH];
    const X_X: PokeFieldI = PokeFieldI::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldI = PokeFieldI::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    const FULL_TWO_TORSION_EXP: usize = 129;
    const EFFECTIVE_TWO_TORSION_EXP: usize = FULL_TWO_TORSION_EXP - 2;
    const THREE_TORSION_EXP: usize = 164;
    const FIVE_TORSION_EXP: usize = 18;

    pub fn get_params() -> PublicParams<PokeFieldI> {
        let cofactor = BigNum::one();

        let effective_two_torsion_order = BigNum::from_prime_power(2, EFFECTIVE_TWO_TORSION_EXP);
        let full_two_torsion_order = 4 * &effective_two_torsion_order;
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );

        let three_torsion_order = BigNum::from_prime_power(3, THREE_TORSION_EXP);
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );

        let five_torsion_order = BigNum::from_prime_power(5, FIVE_TORSION_EXP);
        let five_torsion_cofactor = &full_two_torsion_order * &three_torsion_order * &cofactor;
        let five_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&X_X),
            &PointX::from_x_coord(&Y_X),
            &PointX::from_x_coord(&XY_X),
        );

        let five_adic_basis = (0..=FIVE_TORSION_EXP)
            .map(|exp| BigNum::from_prime_power(5, exp))
            .collect::<Vec<_>>();

        PublicParams {
            starting_curve: Curve::new(&PokeFieldI::from_i32(0)),
            full_two_torsion_exp: FULL_TWO_TORSION_EXP,
            full_two_torsion_order,
            effective_two_torsion_exp: EFFECTIVE_TWO_TORSION_EXP,
            effective_two_torsion_order,
            three_torsion_exp: THREE_TORSION_EXP,
            three_torsion_order,
            five_torsion_exp: FIVE_TORSION_EXP,
            five_torsion_order,
            five_torsion_cofactor,
            two_torsion_basis,
            three_torsion_basis,
            five_torsion_basis,
            five_adic_basis,
        }
    }
}

pub mod poke_iii {
    use num_bigint::BigUint;

    use super::*;
    use crate::{
        fields::{PokeFieldIII, PokeFieldIIIBase},
        poke::PublicParams,
    };

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        125, 35, 156, 116, 70, 116, 20, 147, 135, 197, 97, 128, 14, 20, 217, 232, 110, 5, 5, 148,
        186, 147, 141, 98, 78, 74, 128, 247, 245, 31, 198, 53, 202, 67, 67, 59, 85, 37, 237, 235,
        114, 149, 240, 124, 140, 119, 160, 186, 58, 192, 50, 247, 189, 77, 118, 195, 15, 94, 210,
        115, 21, 149, 27, 215, 165, 188, 73, 48, 147, 202, 151, 87, 34, 120, 227, 191, 170, 86,
        188, 60, 56,
    ];
    const P_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        153, 105, 63, 42, 164, 169, 198, 9, 37, 213, 170, 231, 167, 161, 128, 99, 165, 174, 152,
        70, 192, 160, 74, 253, 130, 156, 142, 103, 243, 207, 138, 35, 35, 240, 178, 123, 99, 111,
        77, 7, 46, 81, 209, 214, 132, 125, 112, 118, 166, 15, 195, 43, 213, 97, 238, 45, 171, 139,
        110, 4, 173, 246, 242, 188, 88, 87, 173, 197, 62, 88, 89, 175, 157, 192, 6, 20, 179, 10,
        171, 250, 172,
    ];
    const Q_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        135, 73, 29, 158, 219, 73, 102, 93, 171, 181, 5, 228, 223, 242, 238, 98, 219, 214, 100,
        131, 154, 205, 157, 50, 166, 236, 202, 148, 191, 24, 44, 173, 158, 133, 26, 110, 18, 98,
        207, 19, 107, 236, 244, 53, 149, 218, 49, 149, 28, 229, 29, 184, 228, 16, 149, 103, 193,
        235, 17, 216, 62, 253, 202, 80, 144, 199, 181, 82, 103, 194, 246, 137, 112, 198, 152, 187,
        240, 135, 181, 188, 24,
    ];
    const Q_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        90, 59, 168, 223, 226, 163, 150, 31, 8, 120, 247, 211, 254, 148, 161, 196, 42, 183, 128,
        121, 92, 218, 156, 204, 149, 160, 19, 21, 171, 163, 78, 229, 184, 165, 5, 120, 107, 200,
        122, 84, 244, 20, 124, 116, 158, 11, 188, 175, 75, 47, 145, 45, 68, 3, 235, 125, 171, 215,
        86, 77, 44, 25, 8, 111, 165, 211, 126, 51, 241, 32, 74, 89, 47, 153, 239, 84, 131, 234,
        255, 135, 207,
    ];
    const PQ_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        226, 209, 57, 13, 100, 243, 93, 90, 76, 108, 96, 91, 209, 84, 16, 207, 48, 76, 114, 135,
        26, 44, 12, 171, 252, 76, 199, 246, 213, 231, 126, 89, 241, 160, 179, 198, 111, 97, 9, 244,
        103, 121, 180, 75, 221, 239, 102, 238, 202, 158, 103, 88, 245, 63, 253, 66, 163, 208, 233,
        118, 105, 242, 255, 76, 133, 71, 158, 101, 101, 72, 88, 216, 135, 128, 240, 236, 59, 107,
        190, 222, 171,
    ];
    const PQ_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        23, 122, 217, 79, 238, 69, 114, 162, 87, 23, 197, 34, 46, 51, 161, 188, 24, 211, 75, 103,
        138, 34, 32, 156, 47, 23, 127, 121, 226, 193, 115, 66, 194, 235, 36, 104, 138, 62, 217,
        186, 122, 167, 123, 9, 107, 155, 248, 182, 89, 106, 134, 131, 117, 94, 207, 243, 52, 227,
        179, 49, 3, 210, 102, 47, 253, 144, 243, 139, 131, 23, 174, 201, 251, 250, 17, 72, 181,
        234, 69, 19, 89,
    ];
    const P_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        139, 128, 56, 106, 219, 186, 155, 228, 57, 84, 54, 21, 151, 154, 80, 58, 175, 50, 43, 97,
        107, 145, 218, 19, 157, 204, 42, 227, 246, 231, 64, 86, 71, 148, 160, 207, 187, 88, 109,
        122, 90, 167, 165, 31, 252, 130, 144, 24, 203, 151, 118, 193, 28, 176, 123, 62, 253, 43,
        120, 234, 216, 194, 236, 254, 106, 230, 18, 73, 1, 218, 84, 199, 106, 96, 102, 155, 48,
        134, 225, 199, 198,
    ];
    const R_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        209, 146, 122, 121, 192, 206, 23, 161, 63, 78, 62, 236, 249, 97, 204, 249, 186, 143, 126,
        210, 232, 46, 57, 92, 17, 222, 124, 140, 29, 137, 254, 144, 61, 145, 129, 0, 198, 95, 95,
        78, 155, 251, 14, 67, 90, 121, 234, 37, 3, 196, 3, 68, 174, 216, 31, 117, 33, 75, 102, 171,
        164, 26, 156, 132, 218, 151, 134, 89, 56, 219, 150, 227, 234, 178, 63, 206, 151, 28, 52,
        105, 43,
    ];
    const S_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        153, 3, 235, 75, 84, 214, 90, 224, 198, 109, 186, 140, 23, 105, 20, 76, 2, 185, 255, 31,
        208, 58, 115, 102, 180, 159, 208, 138, 76, 233, 75, 195, 112, 151, 14, 134, 210, 23, 92,
        114, 10, 86, 133, 99, 219, 26, 20, 113, 199, 114, 133, 36, 112, 241, 141, 60, 72, 130, 118,
        223, 170, 54, 14, 162, 71, 245, 45, 136, 187, 97, 118, 49, 213, 57, 59, 176, 141, 226, 148,
        96, 24,
    ];
    const S_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        116, 214, 179, 242, 52, 28, 114, 226, 130, 150, 56, 174, 114, 72, 88, 139, 225, 85, 70, 33,
        55, 72, 246, 255, 158, 138, 119, 127, 157, 234, 71, 192, 63, 61, 222, 141, 57, 139, 204,
        239, 8, 118, 221, 188, 72, 39, 118, 104, 32, 199, 161, 19, 232, 204, 26, 5, 128, 99, 232,
        184, 49, 40, 109, 151, 22, 91, 249, 115, 78, 117, 146, 253, 232, 166, 114, 10, 1, 78, 94,
        104, 151,
    ];
    const RS_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        231, 141, 220, 45, 172, 109, 53, 151, 117, 191, 134, 167, 3, 72, 42, 87, 35, 96, 185, 34,
        126, 116, 219, 76, 55, 167, 199, 253, 87, 177, 16, 195, 47, 100, 100, 230, 198, 189, 189,
        108, 230, 236, 99, 156, 237, 60, 109, 181, 21, 216, 105, 224, 69, 44, 212, 50, 224, 218,
        57, 54, 254, 36, 80, 218, 169, 147, 18, 248, 185, 43, 4, 241, 39, 38, 196, 22, 15, 80, 62,
        175, 69,
    ];
    const RS_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        156, 23, 90, 133, 35, 240, 167, 204, 93, 81, 85, 174, 252, 113, 27, 142, 127, 71, 4, 212,
        56, 33, 0, 163, 232, 196, 90, 205, 202, 149, 190, 196, 131, 242, 23, 82, 10, 101, 100, 248,
        27, 197, 49, 111, 147, 113, 135, 123, 34, 133, 14, 224, 165, 190, 159, 5, 138, 200, 202,
        220, 239, 222, 168, 71, 178, 113, 8, 83, 165, 42, 155, 95, 171, 5, 173, 37, 159, 104, 146,
        144, 182,
    ];
    const R_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Construct basis points for the 5^c-torsion on E_0
    const X_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        2, 174, 111, 96, 234, 6, 90, 65, 92, 104, 221, 43, 155, 58, 159, 211, 54, 222, 189, 115,
        213, 62, 171, 115, 186, 208, 229, 180, 79, 230, 140, 193, 244, 232, 50, 149, 161, 213, 94,
        132, 226, 180, 127, 186, 231, 248, 167, 235, 98, 216, 228, 206, 147, 247, 181, 235, 218,
        92, 167, 3, 21, 22, 12, 148, 60, 124, 28, 171, 130, 127, 242, 77, 244, 198, 65, 58, 244,
        102, 125, 221, 119,
    ];
    const X_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        93, 205, 196, 61, 237, 193, 17, 87, 174, 180, 26, 23, 147, 126, 27, 107, 82, 34, 250, 121,
        199, 54, 78, 57, 103, 216, 76, 69, 0, 255, 213, 101, 161, 221, 35, 50, 208, 209, 27, 74,
        134, 200, 141, 91, 164, 16, 224, 182, 207, 116, 152, 15, 26, 20, 179, 112, 113, 98, 139,
        145, 7, 26, 47, 241, 213, 94, 92, 116, 135, 79, 10, 145, 36, 233, 15, 113, 162, 48, 223,
        209, 164,
    ];
    const Y_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        45, 25, 141, 108, 13, 188, 71, 36, 165, 55, 80, 139, 187, 251, 177, 241, 84, 37, 21, 183,
        43, 54, 150, 160, 190, 39, 146, 176, 5, 131, 49, 155, 220, 130, 241, 24, 26, 165, 45, 205,
        194, 165, 122, 116, 97, 81, 156, 165, 28, 228, 83, 13, 70, 36, 31, 41, 223, 247, 14, 130,
        248, 64, 21, 144, 161, 227, 130, 95, 9, 231, 68, 60, 121, 222, 216, 252, 29, 143, 233, 215,
        193,
    ];
    const Y_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        90, 253, 92, 19, 5, 164, 212, 93, 89, 80, 149, 206, 8, 172, 126, 25, 249, 136, 181, 166,
        47, 2, 124, 63, 174, 232, 21, 161, 70, 84, 102, 55, 46, 94, 10, 194, 195, 165, 49, 250,
        173, 18, 139, 68, 147, 101, 138, 247, 82, 215, 153, 47, 63, 109, 73, 1, 130, 81, 198, 75,
        3, 30, 46, 89, 227, 173, 34, 129, 116, 217, 207, 103, 82, 181, 104, 192, 32, 76, 191, 172,
        128,
    ];
    const XY_X_RE: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        175, 223, 49, 8, 227, 128, 242, 6, 234, 45, 115, 210, 51, 209, 22, 235, 65, 54, 112, 209,
        134, 126, 188, 111, 150, 165, 214, 43, 243, 215, 230, 166, 137, 202, 123, 198, 101, 13, 71,
        162, 123, 247, 74, 174, 254, 24, 67, 188, 67, 88, 1, 39, 152, 126, 129, 183, 35, 61, 174,
        239, 14, 42, 228, 198, 2, 141, 215, 160, 57, 238, 134, 187, 127, 198, 3, 132, 84, 171, 48,
        199, 107,
    ];
    const XY_X_IM: [u8; PokeFieldIIIBase::ENCODED_LENGTH] = [
        161, 213, 222, 45, 151, 227, 117, 81, 84, 125, 148, 193, 158, 7, 167, 62, 165, 87, 237,
        154, 199, 102, 137, 160, 252, 46, 133, 137, 71, 36, 186, 35, 241, 168, 195, 27, 217, 236,
        10, 49, 40, 185, 45, 136, 195, 152, 14, 210, 209, 218, 98, 82, 101, 84, 142, 95, 111, 42,
        2, 224, 222, 132, 180, 122, 145, 18, 12, 81, 73, 235, 45, 232, 157, 156, 206, 72, 187, 10,
        170, 71, 184,
    ];
    const X_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldIII = PokeFieldIII::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    const FULL_TWO_TORSION_EXP: usize = 192;
    const EFFECTIVE_TWO_TORSION_EXP: usize = FULL_TWO_TORSION_EXP - 2;
    const THREE_TORSION_EXP: usize = 243;
    const FIVE_TORSION_EXP: usize = 28;

    pub fn get_params() -> PublicParams<PokeFieldIII> {
        let cofactor = BigNum::from_prime_power(7, 2);

        let effective_two_torsion_order = BigNum::from_prime_power(2, EFFECTIVE_TWO_TORSION_EXP);
        let full_two_torsion_order = 4 * &effective_two_torsion_order;
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );

        // FIXME: There seems to be a bug in isogeny::utilities::bn::prime_power_to_bn_vartime()
        let three_torsion_order = BigNum::new(
            &BigUint::from(3u8)
                .pow(THREE_TORSION_EXP as u32)
                .to_u64_digits(),
        );
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );

        let five_torsion_order = BigNum::from_prime_power(5, FIVE_TORSION_EXP);
        let five_torsion_cofactor = &full_two_torsion_order * &three_torsion_order * &cofactor;
        let five_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&X_X),
            &PointX::from_x_coord(&Y_X),
            &PointX::from_x_coord(&XY_X),
        );

        let five_adic_basis = (0..=FIVE_TORSION_EXP)
            .map(|exp| BigNum::from_prime_power(5, exp))
            .collect::<Vec<_>>();

        PublicParams {
            starting_curve: Curve::new(&PokeFieldIII::from_i32(0)),
            full_two_torsion_exp: FULL_TWO_TORSION_EXP,
            full_two_torsion_order,
            effective_two_torsion_exp: EFFECTIVE_TWO_TORSION_EXP,
            effective_two_torsion_order,
            three_torsion_exp: THREE_TORSION_EXP,
            three_torsion_order,
            five_torsion_exp: FIVE_TORSION_EXP,
            five_torsion_order,
            five_torsion_cofactor,
            two_torsion_basis,
            three_torsion_basis,
            five_torsion_basis,
            five_adic_basis,
        }
    }
}

pub mod poke_v {
    use super::*;
    use crate::{
        fields::{PokeFieldV, PokeFieldVBase},
        poke::PublicParams,
    };

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        116, 87, 213, 54, 30, 39, 206, 152, 63, 246, 60, 169, 43, 150, 178, 138, 66, 42, 9, 214,
        92, 151, 28, 58, 249, 111, 83, 105, 52, 137, 88, 237, 92, 71, 173, 165, 104, 163, 75, 212,
        23, 193, 137, 229, 0, 212, 96, 140, 216, 83, 224, 164, 56, 20, 204, 173, 54, 89, 253, 158,
        201, 180, 118, 116, 153, 194, 192, 148, 24, 3, 180, 206, 96, 207, 58, 76, 205, 151, 24, 74,
        25, 182, 172, 5, 113, 72, 59, 159, 113, 71, 85, 1, 27, 96, 253, 32, 127, 16, 78, 199, 83,
        214, 128, 168, 112, 30, 133, 0,
    ];
    const P_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        77, 148, 179, 73, 144, 109, 187, 1, 78, 24, 6, 15, 154, 162, 203, 21, 202, 118, 10, 5, 231,
        8, 215, 246, 122, 237, 168, 148, 155, 49, 156, 207, 248, 225, 136, 114, 228, 248, 79, 26,
        149, 141, 221, 187, 118, 45, 202, 254, 97, 39, 0, 1, 156, 212, 77, 95, 227, 133, 166, 5,
        178, 100, 113, 102, 69, 253, 5, 90, 172, 65, 176, 216, 72, 37, 83, 210, 114, 110, 95, 8,
        26, 248, 226, 22, 125, 57, 160, 81, 171, 149, 120, 226, 189, 33, 217, 157, 248, 6, 59, 108,
        147, 152, 31, 153, 220, 97, 67, 32,
    ];
    const Q_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        161, 209, 248, 49, 78, 25, 108, 83, 242, 235, 21, 56, 120, 232, 248, 252, 192, 57, 156,
        140, 127, 26, 137, 194, 68, 103, 254, 154, 226, 214, 212, 218, 8, 223, 54, 169, 193, 40,
        182, 253, 188, 254, 0, 209, 104, 61, 235, 131, 241, 102, 2, 54, 106, 1, 252, 106, 15, 120,
        197, 15, 108, 45, 36, 123, 21, 137, 42, 170, 40, 36, 247, 197, 163, 6, 176, 23, 136, 168,
        25, 56, 84, 111, 175, 36, 119, 41, 24, 225, 209, 238, 53, 208, 194, 77, 90, 21, 56, 237,
        126, 127, 154, 164, 195, 227, 176, 202, 186, 31,
    ];
    const Q_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        111, 224, 111, 218, 159, 193, 114, 125, 247, 51, 231, 28, 34, 98, 198, 67, 33, 220, 60,
        192, 225, 18, 32, 143, 243, 75, 156, 147, 73, 207, 233, 95, 223, 156, 82, 83, 159, 141,
        227, 248, 107, 244, 140, 3, 206, 174, 193, 5, 41, 16, 201, 78, 79, 190, 226, 139, 251, 146,
        218, 10, 130, 120, 104, 189, 241, 178, 98, 158, 3, 85, 205, 71, 80, 57, 22, 170, 152, 210,
        51, 88, 13, 123, 137, 64, 226, 103, 62, 105, 70, 213, 167, 168, 234, 150, 44, 152, 31, 217,
        189, 126, 45, 184, 56, 7, 54, 105, 247, 32,
    ];
    const PQ_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        167, 51, 108, 240, 181, 170, 131, 50, 199, 175, 5, 90, 37, 101, 231, 105, 90, 205, 119, 28,
        57, 225, 87, 202, 224, 141, 189, 71, 133, 27, 72, 250, 10, 17, 143, 126, 2, 122, 222, 64,
        106, 84, 107, 138, 34, 69, 55, 198, 225, 52, 153, 203, 55, 227, 198, 45, 184, 193, 139, 31,
        230, 250, 158, 209, 146, 4, 34, 181, 61, 119, 158, 85, 234, 1, 114, 231, 190, 5, 98, 47,
        11, 182, 165, 113, 36, 19, 172, 0, 221, 27, 127, 20, 244, 240, 84, 209, 59, 114, 121, 97,
        114, 211, 85, 139, 228, 118, 53, 44,
    ];
    const PQ_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        120, 70, 36, 95, 86, 52, 19, 223, 94, 0, 65, 28, 147, 25, 85, 220, 196, 189, 153, 40, 233,
        109, 116, 115, 68, 105, 100, 21, 208, 250, 185, 51, 23, 224, 31, 50, 118, 189, 91, 88, 25,
        84, 156, 123, 191, 106, 49, 192, 55, 232, 209, 82, 18, 222, 130, 85, 4, 1, 209, 95, 116,
        255, 238, 149, 91, 110, 235, 113, 63, 118, 133, 145, 131, 171, 63, 47, 5, 58, 241, 198,
        103, 73, 17, 159, 191, 6, 139, 176, 151, 216, 95, 16, 173, 125, 233, 153, 207, 35, 148,
        134, 44, 25, 127, 5, 150, 119, 109, 42,
    ];
    const P_X: PokeFieldV = PokeFieldV::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldV = PokeFieldV::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldV = PokeFieldV::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        122, 232, 9, 53, 39, 77, 85, 171, 214, 108, 98, 19, 66, 144, 168, 32, 186, 203, 142, 117,
        59, 60, 250, 184, 8, 39, 206, 168, 245, 1, 101, 213, 14, 153, 240, 81, 144, 197, 68, 200,
        155, 37, 206, 94, 134, 179, 127, 101, 99, 251, 53, 60, 78, 5, 110, 214, 43, 211, 63, 103,
        188, 155, 149, 163, 11, 112, 153, 74, 199, 56, 12, 26, 84, 72, 5, 40, 109, 181, 171, 107,
        193, 212, 10, 195, 239, 96, 35, 212, 245, 178, 255, 99, 166, 188, 156, 234, 157, 63, 160,
        38, 106, 155, 178, 255, 0, 93, 111, 70,
    ];
    const R_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        157, 163, 49, 31, 10, 180, 242, 99, 42, 43, 208, 253, 218, 4, 119, 166, 66, 36, 244, 218,
        191, 12, 108, 221, 17, 228, 170, 72, 64, 207, 39, 18, 47, 98, 41, 35, 46, 47, 218, 187,
        228, 250, 125, 233, 94, 220, 184, 252, 247, 164, 209, 158, 100, 176, 48, 27, 207, 117, 95,
        100, 140, 66, 221, 102, 67, 5, 59, 243, 201, 241, 198, 107, 122, 242, 197, 50, 106, 100, 2,
        79, 225, 144, 78, 221, 53, 222, 59, 58, 41, 55, 101, 186, 233, 85, 72, 14, 121, 184, 224,
        205, 214, 167, 21, 91, 84, 178, 3, 18,
    ];
    const S_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        112, 156, 60, 201, 138, 252, 151, 228, 123, 41, 198, 144, 60, 95, 51, 112, 116, 244, 152,
        17, 175, 40, 162, 64, 46, 46, 185, 158, 32, 45, 111, 176, 230, 110, 145, 77, 5, 119, 245,
        173, 136, 70, 77, 11, 231, 243, 159, 175, 219, 69, 5, 18, 163, 207, 122, 57, 244, 226, 122,
        244, 44, 84, 152, 15, 228, 192, 173, 103, 14, 81, 115, 95, 193, 54, 94, 239, 152, 219, 132,
        115, 127, 66, 113, 65, 4, 215, 10, 9, 179, 164, 157, 199, 25, 8, 104, 145, 23, 40, 123,
        185, 170, 99, 13, 10, 157, 14, 6, 57,
    ];
    const S_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        92, 18, 100, 185, 34, 89, 193, 77, 62, 163, 239, 84, 98, 158, 35, 13, 49, 100, 173, 207,
        109, 17, 193, 113, 152, 101, 237, 97, 189, 133, 50, 124, 223, 85, 36, 8, 211, 224, 167,
        199, 5, 104, 213, 200, 251, 34, 252, 117, 94, 143, 72, 135, 204, 21, 49, 168, 213, 80, 165,
        161, 169, 202, 96, 68, 56, 204, 215, 107, 145, 93, 98, 77, 94, 54, 152, 26, 253, 1, 111,
        18, 38, 84, 94, 37, 251, 17, 73, 202, 250, 120, 235, 251, 208, 155, 179, 36, 27, 143, 154,
        200, 82, 15, 10, 177, 41, 153, 107, 59,
    ];
    const RS_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        171, 202, 176, 93, 231, 165, 65, 239, 109, 189, 68, 80, 154, 2, 144, 129, 27, 188, 202, 38,
        67, 41, 58, 168, 118, 120, 8, 232, 204, 19, 78, 247, 86, 101, 251, 245, 204, 36, 139, 22,
        63, 200, 182, 35, 213, 212, 219, 160, 167, 131, 35, 64, 48, 236, 94, 63, 32, 175, 76, 126,
        106, 81, 196, 161, 225, 229, 158, 145, 45, 30, 16, 226, 28, 46, 168, 124, 176, 138, 240,
        243, 2, 47, 68, 17, 106, 164, 188, 65, 77, 28, 114, 122, 33, 222, 177, 243, 31, 149, 240,
        34, 2, 86, 70, 155, 148, 244, 141, 66,
    ];
    const RS_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        49, 187, 217, 101, 242, 197, 87, 135, 141, 30, 73, 52, 32, 195, 10, 224, 41, 8, 20, 88,
        148, 127, 7, 202, 65, 155, 33, 254, 149, 238, 107, 137, 96, 47, 156, 3, 69, 75, 235, 105,
        49, 58, 99, 105, 149, 28, 199, 213, 73, 29, 3, 144, 136, 104, 152, 241, 135, 80, 179, 152,
        78, 115, 129, 11, 188, 149, 202, 235, 69, 61, 167, 88, 145, 250, 111, 34, 234, 236, 19,
        112, 20, 26, 116, 139, 64, 222, 245, 234, 147, 79, 188, 9, 10, 47, 118, 19, 40, 125, 121,
        123, 72, 183, 194, 106, 90, 135, 159, 33,
    ];
    const R_X: PokeFieldV = PokeFieldV::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldV = PokeFieldV::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldV = PokeFieldV::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Construct basis points for the 5^c-torsion on E_0
    const X_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        232, 173, 126, 39, 73, 160, 50, 68, 126, 47, 164, 155, 77, 147, 243, 224, 213, 204, 68, 68,
        99, 221, 143, 1, 18, 14, 47, 230, 67, 160, 45, 4, 162, 62, 141, 41, 212, 176, 137, 69, 120,
        189, 58, 31, 151, 144, 190, 102, 167, 159, 6, 98, 96, 119, 65, 11, 36, 194, 170, 109, 53,
        105, 119, 195, 229, 63, 252, 122, 70, 60, 184, 197, 217, 151, 73, 73, 199, 108, 195, 195,
        228, 8, 11, 233, 183, 48, 122, 176, 231, 186, 81, 119, 160, 252, 224, 15, 56, 108, 77, 183,
        178, 51, 54, 45, 78, 248, 203, 15,
    ];
    const X_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        108, 187, 203, 226, 95, 91, 70, 114, 130, 140, 239, 195, 167, 168, 79, 227, 219, 21, 33,
        50, 231, 211, 31, 206, 49, 6, 150, 6, 186, 124, 25, 248, 246, 215, 189, 175, 200, 57, 211,
        213, 34, 168, 112, 188, 125, 15, 210, 118, 71, 38, 81, 219, 222, 176, 209, 90, 47, 246,
        140, 87, 72, 255, 111, 24, 99, 18, 142, 42, 141, 35, 62, 204, 255, 113, 227, 78, 197, 195,
        251, 137, 217, 161, 18, 211, 214, 6, 109, 129, 244, 160, 235, 18, 182, 177, 187, 208, 229,
        21, 185, 128, 32, 126, 113, 125, 118, 100, 140, 17,
    ];
    const Y_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        165, 205, 168, 250, 252, 98, 207, 168, 106, 173, 107, 145, 50, 222, 24, 131, 81, 185, 128,
        219, 209, 191, 20, 120, 166, 11, 37, 55, 140, 136, 210, 27, 49, 158, 145, 120, 56, 144, 69,
        13, 213, 229, 199, 199, 120, 232, 219, 22, 2, 75, 252, 226, 46, 79, 116, 50, 3, 95, 28,
        108, 181, 161, 105, 63, 3, 4, 124, 100, 235, 212, 43, 205, 26, 79, 163, 60, 187, 137, 101,
        177, 126, 48, 22, 112, 156, 63, 212, 192, 10, 211, 40, 1, 242, 42, 216, 169, 192, 128, 142,
        2, 17, 60, 160, 43, 234, 16, 39, 47,
    ];
    const Y_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        135, 180, 49, 128, 49, 9, 151, 64, 186, 142, 37, 108, 144, 84, 33, 0, 55, 249, 132, 6, 216,
        98, 93, 109, 75, 37, 247, 249, 20, 249, 38, 226, 210, 213, 206, 226, 21, 116, 120, 152,
        175, 74, 26, 238, 77, 29, 141, 118, 157, 0, 54, 209, 175, 14, 3, 98, 253, 172, 201, 240,
        46, 102, 11, 158, 185, 42, 170, 5, 111, 27, 45, 136, 155, 195, 145, 152, 211, 89, 211, 237,
        185, 243, 39, 85, 247, 209, 136, 32, 56, 120, 197, 138, 240, 22, 17, 202, 127, 2, 251, 132,
        26, 47, 94, 18, 100, 12, 15, 51,
    ];
    const XY_X_RE: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        201, 70, 246, 190, 226, 166, 127, 131, 243, 78, 202, 30, 64, 99, 131, 254, 168, 59, 64,
        200, 226, 187, 10, 23, 18, 211, 82, 200, 92, 111, 67, 99, 200, 74, 30, 43, 4, 187, 253,
        233, 53, 193, 62, 242, 68, 159, 176, 88, 205, 194, 206, 159, 42, 156, 208, 174, 154, 117,
        97, 132, 224, 12, 86, 66, 246, 158, 31, 102, 76, 249, 44, 26, 54, 111, 176, 135, 77, 186,
        10, 125, 209, 160, 231, 222, 15, 138, 202, 196, 69, 154, 225, 38, 22, 198, 210, 74, 206,
        117, 121, 117, 209, 152, 5, 99, 190, 46, 201, 47,
    ];
    const XY_X_IM: [u8; PokeFieldVBase::ENCODED_LENGTH] = [
        32, 90, 154, 219, 192, 87, 189, 253, 51, 160, 204, 61, 208, 59, 0, 61, 74, 82, 119, 163, 4,
        158, 148, 178, 102, 150, 101, 236, 181, 11, 208, 228, 29, 189, 173, 210, 53, 242, 166, 58,
        74, 125, 80, 110, 57, 95, 76, 98, 66, 234, 193, 43, 132, 176, 203, 228, 3, 231, 203, 64,
        93, 7, 143, 213, 17, 226, 128, 10, 30, 238, 199, 178, 112, 30, 63, 189, 96, 75, 199, 251,
        138, 82, 92, 35, 231, 186, 156, 35, 143, 176, 2, 48, 4, 45, 197, 186, 0, 69, 36, 61, 103,
        252, 83, 131, 91, 43, 196, 45,
    ];
    const X_X: PokeFieldV = PokeFieldV::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldV = PokeFieldV::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldV = PokeFieldV::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    const FULL_TWO_TORSION_EXP: usize = 256;
    const EFFECTIVE_TWO_TORSION_EXP: usize = FULL_TWO_TORSION_EXP - 2;
    const THREE_TORSION_EXP: usize = 324;
    const FIVE_TORSION_EXP: usize = 36;

    pub fn get_params() -> PublicParams<PokeFieldV> {
        let cofactor = BigNum::from_prime(547);

        let effective_two_torsion_order = BigNum::from_prime_power(2, EFFECTIVE_TWO_TORSION_EXP);
        let full_two_torsion_order = 4 * &effective_two_torsion_order;
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );

        let three_torsion_order = BigNum::from_prime_power(3, THREE_TORSION_EXP);
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );

        let five_torsion_order = BigNum::from_prime_power(5, FIVE_TORSION_EXP);
        let five_torsion_cofactor = &full_two_torsion_order * &three_torsion_order * &cofactor;
        let five_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&X_X),
            &PointX::from_x_coord(&Y_X),
            &PointX::from_x_coord(&XY_X),
        );

        let five_adic_basis = (0..=FIVE_TORSION_EXP)
            .map(|exp| BigNum::from_prime_power(5, exp))
            .collect::<Vec<_>>();

        PublicParams {
            starting_curve: Curve::new(&PokeFieldV::from_i32(0)),
            full_two_torsion_exp: FULL_TWO_TORSION_EXP,
            full_two_torsion_order,
            effective_two_torsion_exp: EFFECTIVE_TWO_TORSION_EXP,
            effective_two_torsion_order,
            three_torsion_exp: THREE_TORSION_EXP,
            three_torsion_order,
            five_torsion_exp: FIVE_TORSION_EXP,
            five_torsion_order,
            five_torsion_cofactor,
            two_torsion_basis,
            three_torsion_basis,
            five_torsion_basis,
            five_adic_basis,
        }
    }
}

pub mod inke_i {
    use super::*;
    use crate::{
        fields::{InkeFieldI, InkeFieldIBase},
        inke::PublicParams,
    };

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        192, 183, 86, 149, 226, 222, 255, 120, 182, 237, 206, 70, 164, 6, 47, 239, 247, 203, 194,
        222, 95, 49, 2, 0, 70, 143, 39, 23, 189, 230, 19, 227, 51, 67, 69, 75, 44, 206, 222, 99,
        117, 77, 90, 186, 16, 104, 191, 147, 98,
    ];
    const P_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        103, 20, 45, 159, 190, 74, 229, 27, 40, 110, 246, 167, 114, 121, 59, 147, 227, 222, 89,
        200, 106, 18, 254, 33, 35, 144, 213, 47, 33, 61, 242, 61, 91, 210, 85, 208, 44, 7, 174,
        197, 57, 51, 57, 248, 222, 45, 22, 51, 51,
    ];
    const Q_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        192, 183, 86, 149, 226, 222, 255, 120, 182, 237, 206, 70, 164, 6, 47, 239, 247, 203, 194,
        222, 95, 49, 2, 0, 70, 143, 39, 23, 189, 230, 19, 227, 51, 67, 69, 75, 44, 206, 222, 99,
        117, 77, 90, 186, 16, 104, 191, 147, 98,
    ];
    const Q_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        152, 235, 210, 96, 65, 181, 26, 228, 215, 145, 9, 88, 141, 134, 196, 108, 19, 233, 136,
        172, 43, 32, 150, 155, 88, 242, 136, 34, 104, 154, 48, 198, 24, 14, 80, 84, 164, 56, 21,
        226, 167, 4, 226, 165, 186, 151, 174, 117, 164,
    ];
    const PQ_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        199, 25, 149, 66, 206, 100, 5, 92, 190, 57, 147, 159, 33, 79, 218, 137, 166, 189, 178, 41,
        212, 25, 119, 242, 153, 215, 95, 201, 165, 211, 241, 142, 171, 154, 39, 20, 141, 171, 37,
        161, 182, 218, 227, 33, 210, 44, 192, 30, 210,
    ];
    const PQ_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [0; InkeFieldIBase::ENCODED_LENGTH];
    const P_X: InkeFieldI = InkeFieldI::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: InkeFieldI = InkeFieldI::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: InkeFieldI = InkeFieldI::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        108, 167, 135, 111, 178, 224, 28, 174, 96, 181, 136, 82, 249, 187, 46, 204, 206, 49, 3,
        204, 114, 40, 146, 39, 5, 71, 82, 132, 158, 215, 28, 134, 196, 214, 95, 88, 230, 233, 133,
        99, 152, 170, 38, 65, 82, 174, 219, 91, 4,
    ];
    const R_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        208, 75, 131, 19, 229, 75, 213, 65, 93, 141, 216, 128, 54, 73, 161, 65, 25, 204, 153, 123,
        0, 157, 106, 204, 5, 203, 106, 228, 90, 94, 218, 91, 75, 82, 44, 119, 210, 239, 91, 174,
        101, 137, 241, 236, 172, 199, 222, 84, 96,
    ];
    const S_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        108, 167, 135, 111, 178, 224, 28, 174, 96, 181, 136, 82, 249, 187, 46, 204, 206, 49, 3,
        204, 114, 40, 146, 39, 5, 71, 82, 132, 158, 215, 28, 134, 196, 214, 95, 88, 230, 233, 133,
        99, 152, 170, 38, 65, 82, 174, 219, 91, 4,
    ];
    const S_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        47, 180, 124, 236, 26, 180, 42, 190, 162, 114, 39, 127, 201, 182, 94, 190, 221, 251, 72,
        249, 149, 149, 41, 241, 117, 183, 243, 109, 46, 121, 72, 168, 40, 142, 121, 173, 254, 79,
        103, 249, 123, 174, 41, 177, 236, 253, 229, 83, 119,
    ];
    const RS_X_RE: [u8; InkeFieldIBase::ENCODED_LENGTH] = [
        77, 147, 144, 129, 1, 122, 226, 93, 70, 80, 179, 245, 103, 143, 211, 87, 29, 115, 84, 244,
        63, 170, 158, 138, 40, 113, 38, 113, 230, 251, 245, 173, 129, 20, 114, 96, 158, 182, 208,
        161, 75, 71, 228, 191, 250, 107, 148, 107, 19,
    ];
    const RS_X_IM: [u8; InkeFieldIBase::ENCODED_LENGTH] = [0; InkeFieldIBase::ENCODED_LENGTH];
    const R_X: InkeFieldI = InkeFieldI::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: InkeFieldI = InkeFieldI::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: InkeFieldI = InkeFieldI::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    const FULL_TWO_TORSION_EXP: usize = 128;
    const EFFECTIVE_TWO_TORSION_EXP: usize = FULL_TWO_TORSION_EXP - 2;
    const THREE_TORSION_EXP: usize = 162;

    pub fn get_params() -> PublicParams<InkeFieldI> {
        let effective_two_torsion_order = BigNum::from_prime_power(2, EFFECTIVE_TWO_TORSION_EXP);
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );

        let three_torsion_order = BigNum::from_prime_power(3, THREE_TORSION_EXP);
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );

        PublicParams {
            starting_curve: Curve::new(&InkeFieldI::from_i32(0)),
            effective_two_torsion_exp: EFFECTIVE_TWO_TORSION_EXP,
            effective_two_torsion_order,
            three_torsion_exp: THREE_TORSION_EXP,
            three_torsion_order,
            two_torsion_basis,
            three_torsion_basis,
        }
    }
}
