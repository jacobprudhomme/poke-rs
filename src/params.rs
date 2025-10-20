use fp2::traits::Fp as _;
use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX};

use crate::PublicParams;

pub mod poke_i {
    use super::*;
    use crate::fields::PokeFieldI;

    // Construct basis points for the 2^a-torsion on E_0
    const P_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const P_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        156, 1, 55, 113, 211, 191, 79, 224, 97, 54, 107, 37, 254, 167, 210, 138, 199, 125, 108,
        159, 62, 27, 61, 12, 176, 93, 127, 206, 236, 40, 77, 235, 18, 81, 163, 191, 61, 216, 30,
        105, 141, 244, 112, 38, 122, 199, 207, 251, 158, 170, 70, 187, 238, 73,
    ];
    const Q_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        225, 204, 48, 233, 181, 190, 18, 222, 247, 38, 59, 93, 252, 209, 65, 62, 195, 253, 222, 58,
        179, 18, 119, 130, 98, 196, 148, 139, 59, 204, 93, 73, 22, 7, 63, 63, 184, 164, 108, 255,
        205, 79, 133, 20, 182, 27, 46, 205, 220, 82, 131, 215, 39, 28,
    ];
    const Q_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        99, 254, 200, 142, 44, 64, 176, 31, 158, 201, 148, 218, 1, 88, 45, 117, 138, 22, 138, 198,
        255, 79, 247, 48, 32, 60, 177, 74, 18, 196, 166, 147, 100, 195, 53, 43, 189, 187, 224, 237,
        23, 35, 156, 158, 249, 126, 66, 202, 241, 91, 209, 56, 29, 32,
    ];
    const PQ_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        84, 5, 39, 47, 28, 26, 72, 53, 222, 25, 244, 169, 18, 138, 123, 250, 227, 91, 135, 191,
        182, 168, 208, 156, 231, 66, 10, 171, 57, 90, 207, 9, 222, 195, 240, 102, 7, 222, 148, 122,
        208, 175, 249, 130, 55, 245, 12, 92, 175, 174, 252, 231, 208, 15,
    ];
    const PQ_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const P_X: PokeFieldI = PokeFieldI::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldI = PokeFieldI::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Construct basis points for the 3^b-torsion on E_0
    const R_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const R_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        242, 78, 54, 243, 244, 158, 40, 209, 213, 36, 144, 132, 126, 115, 146, 252, 115, 95, 79,
        49, 121, 90, 228, 120, 114, 82, 233, 129, 214, 22, 113, 22, 116, 81, 115, 222, 238, 180,
        157, 29, 159, 205, 134, 216, 253, 65, 214, 79, 148, 149, 147, 24, 195, 7,
    ];
    const S_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        33, 100, 192, 219, 131, 122, 237, 66, 95, 99, 60, 177, 230, 250, 51, 190, 104, 113, 44,
        242, 139, 87, 147, 181, 249, 53, 197, 220, 252, 127, 88, 234, 23, 241, 221, 97, 160, 52,
        102, 44, 37, 165, 139, 203, 245, 120, 204, 216, 248, 102, 186, 121, 47, 14,
    ];
    const S_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        13, 177, 201, 12, 11, 97, 215, 46, 42, 219, 111, 123, 129, 140, 109, 3, 222, 52, 167, 52,
        197, 16, 80, 196, 93, 71, 71, 151, 40, 214, 130, 104, 3, 195, 101, 12, 12, 223, 97, 57, 6,
        74, 134, 236, 117, 4, 60, 118, 252, 112, 132, 219, 72, 98,
    ];
    const RS_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        214, 198, 30, 107, 72, 179, 73, 69, 35, 50, 116, 38, 27, 143, 85, 161, 55, 176, 109, 176,
        64, 247, 227, 127, 52, 115, 253, 72, 217, 177, 78, 213, 224, 192, 75, 192, 253, 45, 130,
        177, 170, 220, 184, 89, 185, 137, 120, 89, 231, 163, 80, 255, 92, 95,
    ];
    const RS_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const R_X: PokeFieldI = PokeFieldI::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldI = PokeFieldI::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldI = PokeFieldI::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Construct basis points for the 5^c-torsion on E_0
    const X_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const X_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        211, 24, 39, 150, 242, 9, 42, 168, 202, 224, 82, 61, 102, 182, 231, 124, 213, 107, 144, 72,
        84, 15, 181, 210, 65, 156, 234, 60, 141, 56, 253, 222, 254, 41, 3, 136, 237, 101, 182, 89,
        189, 117, 17, 158, 8, 209, 192, 197, 185, 255, 80, 133, 107, 67,
    ];
    const Y_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        142, 141, 154, 162, 86, 251, 208, 110, 83, 81, 167, 239, 99, 27, 248, 99, 176, 209, 50, 79,
        95, 226, 187, 103, 115, 94, 168, 239, 128, 125, 222, 127, 12, 58, 148, 85, 96, 16, 38, 236,
        30, 216, 153, 163, 196, 201, 222, 27, 117, 237, 189, 56, 217, 95,
    ];
    const Y_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        44, 231, 216, 105, 13, 246, 213, 87, 53, 31, 173, 194, 153, 73, 24, 131, 124, 40, 102, 29,
        234, 91, 127, 106, 142, 253, 69, 220, 113, 180, 246, 159, 120, 234, 213, 98, 13, 46, 73,
        253, 231, 161, 251, 38, 107, 117, 81, 0, 215, 6, 199, 110, 160, 38,
    ];
    const XY_X_RE: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [
        7, 127, 166, 170, 67, 68, 236, 216, 13, 61, 70, 222, 190, 115, 147, 244, 207, 140, 116,
        141, 195, 61, 63, 202, 239, 236, 93, 15, 92, 242, 111, 151, 53, 67, 144, 196, 218, 77, 91,
        160, 138, 225, 199, 32, 138, 40, 48, 229, 231, 107, 234, 178, 109, 61,
    ];
    const XY_X_IM: [u8; PokeFieldI::ENCODED_LENGTH / 2] = [0; PokeFieldI::ENCODED_LENGTH / 2];
    const X_X: PokeFieldI = PokeFieldI::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldI = PokeFieldI::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    pub fn create_poke_i_params() -> PublicParams<PokeFieldI> {
        let two_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&P_X),
            &PointX::from_x_coord(&Q_X),
            &PointX::from_x_coord(&PQ_X),
        );
        let three_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&R_X),
            &PointX::from_x_coord(&S_X),
            &PointX::from_x_coord(&RS_X),
        );
        let five_torsion_basis = BasisX::from_points(
            &PointX::from_x_coord(&X_X),
            &PointX::from_x_coord(&Y_X),
            &PointX::from_x_coord(&XY_X),
        );

        PublicParams {
            starting_curve: Curve::new(&PokeFieldI::from_i32(6)),
            two_torsion_exp: 129,
            three_torsion_exp: 164,
            five_torsion_exp: 18,
            two_torsion_basis,
            three_torsion_basis,
            five_torsion_basis,
        }
    }
}
