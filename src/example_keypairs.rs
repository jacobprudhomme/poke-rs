use std::marker::PhantomData;

use isogeny::elliptic::{basis::BasisX, curve::Curve, point::PointX};
use num_bigint::BigUint;

use crate::{PrvKey, PubKey, SUCCESS_RETVAL};

pub mod poke_i {
    use crate::fields::{PokeFieldI, PokeFieldIBase};

    use super::*;

    // Parameter for codomain curve E_A (obtained by using POKE keygen from INKE project)
    const A_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        246, 65, 240, 120, 17, 123, 195, 253, 248, 157, 66, 168, 182, 103, 3, 123, 180, 239, 78, 7,
        35, 184, 146, 197, 133, 224, 243, 37, 170, 48, 76, 210, 125, 37, 2, 6, 230, 66, 21, 147,
        34, 54, 248, 17, 51, 172, 179, 84, 54, 155, 219, 123, 83, 53,
    ];
    const A_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        36, 96, 100, 30, 231, 141, 157, 203, 171, 31, 14, 146, 71, 177, 78, 48, 206, 164, 91, 91,
        4, 21, 54, 87, 163, 224, 242, 244, 42, 43, 247, 71, 222, 80, 43, 68, 0, 216, 30, 189, 129,
        209, 76, 84, 207, 41, 83, 72, 62, 131, 217, 18, 128, 22,
    ];
    const A: PokeFieldI = PokeFieldI::const_decode_no_check(&A_RE, &A_IM);

    // Masked images of 2^a-torsion basis for E_0 on E_A (obtained by using POKE keygen from INKE project)
    const P_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        65, 190, 116, 164, 250, 75, 184, 179, 92, 22, 104, 27, 249, 242, 37, 193, 212, 238, 19,
        238, 111, 113, 66, 179, 180, 52, 156, 240, 204, 90, 29, 240, 33, 9, 125, 169, 11, 128, 195,
        130, 38, 37, 21, 159, 241, 1, 62, 144, 178, 103, 86, 44, 11, 48,
    ];
    const P_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        123, 161, 202, 150, 156, 29, 217, 106, 145, 222, 116, 218, 169, 73, 215, 60, 40, 76, 93, 0,
        115, 207, 139, 190, 110, 139, 118, 90, 158, 37, 12, 176, 128, 55, 133, 143, 157, 134, 2,
        134, 9, 139, 216, 92, 243, 32, 137, 62, 64, 243, 20, 239, 176, 31,
    ];
    const Q_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        209, 71, 122, 128, 148, 224, 241, 255, 224, 210, 157, 146, 202, 170, 153, 85, 122, 154,
        186, 197, 211, 112, 0, 240, 63, 92, 234, 224, 96, 185, 46, 177, 244, 129, 49, 59, 173, 193,
        108, 70, 193, 32, 180, 211, 30, 117, 151, 141, 125, 98, 61, 92, 148, 76,
    ];
    const Q_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        211, 115, 156, 59, 132, 75, 255, 229, 82, 107, 73, 121, 148, 170, 171, 1, 15, 188, 83, 60,
        155, 109, 241, 18, 155, 190, 65, 165, 5, 174, 171, 196, 24, 5, 150, 181, 193, 244, 122,
        107, 132, 251, 142, 173, 209, 252, 132, 211, 215, 178, 122, 205, 127, 50,
    ];
    const PQ_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        41, 65, 54, 10, 68, 240, 40, 128, 78, 159, 6, 147, 56, 50, 170, 182, 149, 238, 143, 10, 28,
        27, 28, 14, 198, 188, 95, 93, 30, 245, 219, 92, 226, 51, 76, 82, 70, 92, 39, 219, 245, 31,
        123, 58, 60, 199, 218, 187, 225, 60, 227, 162, 41, 15,
    ];
    const PQ_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        246, 144, 13, 139, 196, 188, 231, 60, 159, 189, 12, 191, 32, 70, 186, 122, 20, 37, 81, 239,
        152, 231, 193, 231, 251, 191, 163, 185, 119, 224, 193, 33, 212, 170, 74, 39, 69, 16, 195,
        123, 197, 195, 169, 82, 62, 166, 44, 81, 37, 50, 44, 202, 66, 8,
    ];
    const P_X: PokeFieldI = PokeFieldI::const_decode_no_check(&P_X_RE, &P_X_IM);
    const Q_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Q_X_RE, &Q_X_IM);
    const PQ_X: PokeFieldI = PokeFieldI::const_decode_no_check(&PQ_X_RE, &PQ_X_IM);

    // Masked images of 3^b-torsion basis for E_0 on E_A (obtained by using POKE keygen from INKE project)
    const R_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        233, 87, 81, 241, 124, 190, 179, 48, 37, 176, 145, 140, 48, 91, 48, 85, 133, 40, 137, 146,
        81, 77, 191, 42, 199, 168, 170, 17, 187, 255, 41, 39, 179, 155, 20, 8, 170, 77, 146, 55,
        70, 9, 226, 158, 44, 253, 143, 95, 182, 35, 183, 209, 121, 18,
    ];
    const R_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        70, 192, 74, 131, 222, 143, 1, 76, 45, 243, 89, 146, 197, 215, 30, 241, 144, 87, 45, 251,
        122, 202, 44, 78, 196, 244, 58, 247, 252, 191, 44, 183, 218, 202, 12, 150, 251, 95, 43, 15,
        149, 9, 109, 67, 133, 206, 209, 112, 141, 31, 77, 103, 3, 64,
    ];
    const S_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        68, 242, 225, 144, 154, 199, 142, 113, 162, 148, 93, 152, 107, 70, 255, 74, 133, 197, 112,
        111, 148, 218, 205, 136, 48, 253, 10, 80, 20, 18, 43, 247, 87, 81, 54, 144, 18, 195, 185,
        217, 107, 251, 189, 19, 28, 194, 187, 90, 83, 26, 150, 157, 11, 80,
    ];
    const S_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        194, 160, 58, 42, 125, 57, 227, 64, 202, 221, 83, 63, 213, 46, 46, 119, 189, 219, 17, 234,
        112, 123, 209, 178, 73, 220, 122, 247, 253, 118, 182, 210, 122, 205, 55, 146, 221, 92, 72,
        250, 152, 102, 90, 249, 237, 215, 167, 244, 2, 69, 38, 87, 141, 49,
    ];
    const RS_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        106, 55, 133, 42, 139, 31, 188, 105, 184, 10, 221, 182, 251, 208, 27, 93, 11, 176, 168, 79,
        190, 70, 48, 37, 150, 34, 155, 254, 78, 80, 141, 170, 92, 135, 130, 76, 82, 76, 224, 181,
        51, 110, 134, 38, 10, 29, 54, 77, 94, 88, 160, 255, 13, 105,
    ];
    const RS_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        98, 93, 67, 98, 47, 27, 129, 92, 187, 137, 1, 102, 74, 131, 158, 233, 106, 45, 22, 91, 167,
        29, 128, 222, 196, 75, 154, 88, 245, 76, 210, 97, 82, 209, 163, 255, 139, 129, 102, 120,
        204, 238, 105, 86, 144, 128, 117, 30, 15, 207, 213, 235, 213, 30,
    ];
    const R_X: PokeFieldI = PokeFieldI::const_decode_no_check(&R_X_RE, &R_X_IM);
    const S_X: PokeFieldI = PokeFieldI::const_decode_no_check(&S_X_RE, &S_X_IM);
    const RS_X: PokeFieldI = PokeFieldI::const_decode_no_check(&RS_X_RE, &RS_X_IM);

    // Masked images of 5^c-torsion basis for E_0 on E_A (obtained by using POKE keygen from INKE project)
    const X_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        169, 7, 118, 3, 100, 8, 202, 141, 152, 53, 82, 170, 134, 86, 193, 123, 34, 155, 123, 205,
        118, 210, 11, 98, 197, 223, 220, 219, 97, 160, 142, 249, 143, 19, 62, 218, 82, 231, 193,
        244, 89, 203, 24, 157, 196, 177, 247, 135, 177, 39, 163, 192, 193, 76,
    ];
    const X_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        186, 123, 7, 183, 80, 147, 126, 25, 55, 222, 196, 154, 20, 136, 124, 202, 89, 254, 238,
        104, 238, 211, 207, 149, 28, 173, 37, 4, 75, 232, 37, 165, 55, 133, 149, 163, 205, 242,
        103, 209, 42, 221, 234, 158, 201, 84, 203, 166, 37, 122, 131, 109, 236, 94,
    ];
    const Y_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        139, 0, 58, 231, 105, 62, 247, 117, 138, 55, 184, 39, 36, 196, 23, 245, 110, 53, 64, 244,
        234, 105, 6, 90, 250, 113, 59, 205, 223, 101, 141, 242, 225, 122, 55, 129, 115, 140, 180,
        155, 238, 127, 143, 31, 50, 113, 146, 63, 61, 20, 189, 63, 100, 69,
    ];
    const Y_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        107, 230, 130, 93, 152, 103, 202, 14, 196, 180, 80, 185, 61, 255, 37, 222, 26, 80, 179,
        156, 174, 230, 101, 198, 238, 50, 9, 99, 148, 115, 225, 117, 254, 190, 201, 44, 87, 103,
        220, 249, 109, 80, 138, 34, 211, 244, 6, 22, 150, 239, 186, 60, 6, 101,
    ];
    const XY_X_RE: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        44, 98, 73, 70, 140, 149, 226, 113, 185, 109, 138, 27, 242, 53, 23, 111, 174, 112, 114, 8,
        249, 1, 210, 106, 76, 72, 25, 42, 251, 110, 200, 176, 84, 137, 231, 237, 146, 221, 153,
        155, 207, 33, 123, 172, 1, 233, 34, 68, 58, 237, 222, 222, 91, 73,
    ];
    const XY_X_IM: [u8; PokeFieldIBase::ENCODED_LENGTH] = [
        233, 136, 217, 148, 49, 33, 199, 108, 142, 243, 105, 143, 146, 224, 240, 74, 82, 192, 60,
        89, 239, 248, 123, 128, 169, 243, 205, 157, 94, 148, 127, 249, 171, 202, 4, 96, 92, 163,
        115, 199, 177, 168, 198, 232, 30, 128, 181, 44, 59, 60, 65, 29, 174, 51,
    ];
    const X_X: PokeFieldI = PokeFieldI::const_decode_no_check(&X_X_RE, &X_X_IM);
    const Y_X: PokeFieldI = PokeFieldI::const_decode_no_check(&Y_X_RE, &Y_X_IM);
    const XY_X: PokeFieldI = PokeFieldI::const_decode_no_check(&XY_X_RE, &XY_X_IM);

    // Degree and scalars that make up private key
    const Q: [u8; 16] = [
        123, 49, 95, 50, 21, 90, 35, 134, 150, 96, 179, 53, 178, 88, 46, 24,
    ];
    const ALPHA: [u8; 16] = [
        87, 203, 66, 73, 11, 220, 246, 170, 80, 88, 83, 143, 197, 146, 95, 22,
    ];
    const BETA: [u8; 16] = [
        73, 192, 105, 175, 242, 125, 166, 222, 106, 181, 37, 51, 66, 87, 226, 7,
    ];
    const DELTA: [u8; 6] = [199, 215, 172, 138, 84, 1];

    pub fn get_pub_key() -> PubKey<PokeFieldI> {
        let E_A = Curve::new(&A);

        let xP = PointX::from_x_coord(&P_X);
        let xQ = PointX::from_x_coord(&Q_X);
        let xPQ = PointX::from_x_coord(&PQ_X);
        assert_eq!(E_A.is_on_curve(&xP.x()), SUCCESS_RETVAL, "P is not on E_A");
        assert_eq!(E_A.is_on_curve(&xQ.x()), SUCCESS_RETVAL, "Q is not on E_A");
        assert_eq!(
            E_A.is_on_curve(&xPQ.x()),
            SUCCESS_RETVAL,
            "P - Q is not on E_A"
        );
        // TODO: check that points not only have given order, but that they also don't have smaller order

        let xR = PointX::from_x_coord(&R_X);
        let xS = PointX::from_x_coord(&S_X);
        let xRS = PointX::from_x_coord(&RS_X);
        assert_eq!(E_A.is_on_curve(&xR.x()), SUCCESS_RETVAL, "R is not on E_A");
        assert_eq!(E_A.is_on_curve(&xS.x()), SUCCESS_RETVAL, "S is not on E_A");
        assert_eq!(
            E_A.is_on_curve(&xRS.x()),
            SUCCESS_RETVAL,
            "R - S is not on E_A"
        );

        let xX = PointX::from_x_coord(&X_X);
        let xY = PointX::from_x_coord(&Y_X);
        let xXY = PointX::from_x_coord(&XY_X);
        assert_eq!(E_A.is_on_curve(&xX.x()), SUCCESS_RETVAL, "X is not on E_A");
        assert_eq!(E_A.is_on_curve(&xY.x()), SUCCESS_RETVAL, "Y is not on E_A");
        assert_eq!(
            E_A.is_on_curve(&xXY.x()),
            SUCCESS_RETVAL,
            "X - Y is not on E_A"
        );

        // Construct public key from above values (E_A, (P_A, Q_A), (R_A, S_A), (X_A, Y_A))
        PubKey {
            codomain_curve: E_A,
            masked_two_torsion_basis_img: BasisX::from_points(&xP, &xQ, &xPQ),
            masked_three_torsion_basis_img: BasisX::from_points(&xR, &xS, &xRS),
            masked_five_torsion_basis_img: BasisX::from_points(&xX, &xY, &xXY),
        }
    }

    pub fn get_prv_key() -> PrvKey<PokeFieldI> {
        PrvKey {
            q: BigUint::from_bytes_le(&Q),
            alpha: BigUint::from_bytes_le(&ALPHA),
            beta: BigUint::from_bytes_le(&BETA),
            delta: BigUint::from_bytes_le(&DELTA),
            _field: PhantomData,
        }
    }
}
