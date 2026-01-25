//! Functions for transforming quaternions to endomorphisms and associated helper functions.
//!
//! Adapted from the KLaPoTi-Rust library: https://github.com/isogeny-klapoti/klapoti-rust.

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{curve::Curve, projective_point::Point};
use rug::{
    Integer,
    integer::{IsPrime, Order},
    rand::RandState,
};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

// Primes up to 101 congruent to 1 (mod 4). Taken from "Ready to SQI?"
const SMALL_PRIMES: [u32; 13] = [2, 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101];
// Primes up to 101 congruent to 3 (mod 4). Taken from "Ready to SQI?"
const SMALL_BAD_PRIMES: [u32; 13] = [3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83];

struct SignedInteger {
    n: Integer,
    is_negative: u32,
}

impl SignedInteger {
    fn zero() -> Self {
        Self {
            n: Integer::ZERO,
            is_negative: 0,
        }
    }

    fn new(n: Integer) -> Self {
        Self {
            n: n.clone(),
            is_negative: if n.signum() == -1 {
                SUCCESS_RETVAL
            } else {
                FAILURE_RETVAL
            },
        }
    }
}

pub struct Quaternion {
    x: SignedInteger,
    y: SignedInteger,
    z: SignedInteger,
    t: SignedInteger,
}

impl Quaternion {
    fn zero() -> Self {
        Self {
            x: SignedInteger::zero(),
            y: SignedInteger::zero(),
            z: SignedInteger::zero(),
            t: SignedInteger::zero(),
        }
    }
}

fn complex_mul(
    (x_re, x_im): (&Integer, &Integer),
    (y_re, y_im): (&Integer, &Integer),
) -> (Integer, Integer) {
    let re = x_re.clone() * y_re - x_im * y_im;
    let im = x_re.clone() * y_im + y_re * x_im;

    (re, im)
}

fn complex_mul_by_complex_power(
    (x_re, x_im): (&Integer, &Integer),
    (y_re, y_im): (&Integer, &Integer),
    exp: u32,
) -> (Integer, Integer) {
    let mut re = Integer::ONE.clone();
    let mut im = Integer::ZERO;

    for i in 1..=u32::BITS {
        (re, im) = complex_mul((&re, &im), (&re, &im));
        if (exp >> (u32::BITS - i)) & 1 == 1 {
            (re, im) = complex_mul((&re, &im), (&y_re, &y_im));
        }
    }
    complex_mul((&re, &im), (&x_re, &x_im))
}

// WARN: p is prime
fn tonelli_shanks(a: &Integer, p: &Integer) -> Integer {
    let TWO = Integer::from(2);

    let mut q = p.clone() - Integer::ONE;

    // SAFETY: Since p is prime, q >= 1 so it will have a 1 bit somewhere
    let exp1 = q.find_one(0).unwrap();
    q >>= exp1;

    // Find generator, a quadratic non-residue
    let mut qnr = TWO.clone();
    while qnr.legendre(&p) != -1 {
        qnr += Integer::ONE;
    }
    // SAFETY: we know the exponent is positive
    let mut z = qnr.pow_mod(&q, &p).unwrap();

    // SAFETY: we know the exponent is positive
    let mut y = a.clone().pow_mod(&q, &p).unwrap();

    let exp2 = (q + 1) >> 1;
    // SAFETY: we know the exponent is positive
    let mut x = a.clone().pow_mod(&exp2, &p).unwrap();
    let mut exp3 = Integer::ONE.clone() << (exp1 - 2);

    for _ in 0..exp1 {
        // SAFETY: we know the exponent is positive
        let b = y.clone().pow_mod(&exp3, &p).unwrap();

        if b == p.clone() - Integer::ONE {
            x *= z.clone();
            x %= p;

            y *= z.clone().square();
            y %= y.clone();
        }

        // SAFETY: we know the exponent is positive
        z = z.pow_mod(&TWO, &p).unwrap();
        exp3 >>= 1;
    }

    x
}

fn sqrt_mod_p(a: &Integer, p: &Integer) -> (Integer, u32) {
    let result: Integer;
    let mut retval = SUCCESS_RETVAL;

    let reduced_a = a.clone() % p;
    if reduced_a.jacobi(p) != 1 {
        retval = FAILURE_RETVAL
    };

    if p.clone() % 4 == 3 {
        let exp = (p.clone() + 1) >> 2;
        // SAFETY: we know exp is positive, so this will never fail
        result = reduced_a.pow_mod(&exp, &p).unwrap();
    } else if p.clone() % 8 == 5 {
        let mut exp = (p.clone() - 1) >> 2;
        // SAFETY: we know exp is positive, so this will never fail
        let mut tmp = reduced_a.clone().pow_mod(&exp, &p).unwrap();

        if &tmp == Integer::ONE {
            exp = (p.clone() + 3) >> 3;
            // SAFETY: we know exp is positive, so this will never fail
            result = reduced_a.pow_mod(&exp, &p).unwrap();
        } else {
            exp = (p.clone() - 5) >> 3;

            tmp = 4 * reduced_a.clone();
            // SAFETY: we know exp is positive, so this will never fail
            tmp = tmp.pow_mod(&exp, &p).unwrap();
            tmp = 2 * reduced_a * tmp;

            result = tmp % p;
        }
    } else {
        result = tonelli_shanks(a, p);
    }

    (result, retval)
}

/// Compute `x` and `y` such that x^2 + n*y^2 = p where p is prime.
/// It assumes that there is a sqrt of -1 mod p.
pub fn cornacchia(n: &Integer, p: &Integer) -> (Integer, Integer, u32) {
    let mut retval = SUCCESS_RETVAL;

    let TWO = Integer::from(2);

    if p == &TWO {
        if n == Integer::ONE {
            return (Integer::ONE.clone(), Integer::ONE.clone(), SUCCESS_RETVAL);
        } else {
            return (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL);
        }
    }

    if n.clone() % p == Integer::ZERO {
        return (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL);
    }

    let (mut r2, ok) = sqrt_mod_p(&(-n.clone()), &p);
    retval &= ok;

    let mut r0 = Integer::ZERO;
    let mut r1 = p.clone();
    let mut product = p.clone();

    while &product >= p {
        (_, r0) = r2.div_rem(r1.clone());
        r2 = r1.clone();
        r1 = r0.clone();
        product = r0.clone().square();
    }

    // Test if result is solution
    let mut a = p - product.clone();
    (a, r2) = a.div_rem(n.clone());

    let x: Integer;
    let y = a.sqrt();
    if &r2 == &Integer::ZERO {
        x = r0.clone();
        a = y.clone().square();
        a *= n;
        product = product + a;

        if &product == p {
            return (x.abs(), y.abs(), retval);
        }
    }

    (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL)
}

fn cornacchia_extended_prime_loop(
    c: (&Integer, &Integer),
    p: &Integer,
    valuation_of_p: u32,
) -> (Integer, Integer, u32) {
    let mut retval = SUCCESS_RETVAL;

    let (re, im, ok) = cornacchia(Integer::ONE, p);
    retval &= ok;

    if re != 0 && im != 0 {
        let (re, im) = complex_mul_by_complex_power(c, (&re, &im), valuation_of_p);
        (re, im, retval)
    } else {
        (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL)
    }
}

fn cornacchia_extended(
    n: &Integer,
    primes: &[Integer],
    bad_primes: &[Integer],
) -> (Integer, Integer, u32) {
    let mut retval = SUCCESS_RETVAL;

    if bad_primes
        .iter()
        .any(|prime| n.clone() % prime == Integer::ZERO)
    {
        return (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL);
    }

    let mut n = n.to_owned();
    let mut valuations = vec![0; primes.len()];
    for (i, (p, valuation_at_p)) in primes.iter().zip(valuations.iter_mut()).enumerate() {
        if (p.clone() % 4) == 1 || i == 0 {
            *valuation_at_p = n.remove_factor_mut(p);
        }
    }

    if n.clone() % 4 == 1 {
        let primality_of_remaining_n = n.is_probably_prime(30);
        if primality_of_remaining_n == IsPrime::Probably
            || primality_of_remaining_n == IsPrime::Yes
            || n == 1
        {
            let (mut x, mut y) = if n == 1 {
                (Integer::ONE.clone(), Integer::ZERO)
            } else {
                let (x, y, ok) = cornacchia(Integer::ONE, &n);
                retval &= ok;
                (x, y)
            };

            if !(x == 0 && y == 0) {
                for (p, valuation_at_p) in primes.iter().zip(valuations) {
                    if valuation_at_p != 0 {
                        let (x_, y_, ok) =
                            cornacchia_extended_prime_loop((&x, &y), p, valuation_at_p);
                        (x, y) = (x_, y_);
                        retval &= ok;
                    }
                }
                return (x, y, retval);
            }
        }
    }

    (Integer::ZERO, Integer::ZERO, FAILURE_RETVAL)
}

// Finds a quaternion of given norm in O_0
pub fn represent_integer<const NUM_WORDS: usize, const NUM_WORDS_P: usize>(
    target_norm: &BigNum<NUM_WORDS>,
    p: &BigNum<NUM_WORDS_P>,
) -> (Quaternion, u32) {
    let mut rng = RandState::new();

    let TWO = Integer::from(2);
    let PRIMES = SMALL_PRIMES.map(|prime| Integer::from(prime));
    let BAD_PRIMES = SMALL_BAD_PRIMES.map(|prime| Integer::from(prime));

    let target_norm = Integer::from_digits(&target_norm.to_le_bytes(), Order::Lsf);
    let p = Integer::from_digits(&p.to_le_bytes(), Order::Lsf);
    let target_norm_over_p = &target_norm / p.clone();

    // FIXME: does this do rounding division? And does rounding division mess up our parameters?
    let m_prime = target_norm_over_p.clone().sqrt();
    let m_prime_range = &TWO * m_prime.clone() + Integer::ONE;

    let mut x = Integer::ZERO;
    let mut y = Integer::ZERO;
    let mut z = Integer::ZERO;
    let mut t = Integer::ZERO;

    let mut counter = 1_000_000;
    let mut found_solution = FAILURE_RETVAL;
    while counter > 0 && found_solution == FAILURE_RETVAL {
        z = m_prime_range.clone().random_below(&mut rng) - &m_prime;

        let m_dblprime = (&target_norm_over_p - z.clone().square()).sqrt();
        let m_dblprime_range = &TWO * m_dblprime.clone() + Integer::ONE;
        t = m_dblprime_range.clone().random_below(&mut rng) - &m_dblprime;

        let norm = &target_norm - &p * (z.clone().square() + t.clone().square());
        (x, y, found_solution) = cornacchia_extended(&norm, &PRIMES, &BAD_PRIMES);

        counter -= 1;
    }
    if counter == 0 {
        println!("Tried 1,000,000 times to solve with Cornacchia, but couldn't find a solution");
        return (Quaternion::zero(), FAILURE_RETVAL);
    }

    // Quaternion (x + y*i + z*j + t*k) defined over the basis (1,i,j,k) = (1,i,j,ij)
    (
        Quaternion {
            x: SignedInteger::new(x),
            y: SignedInteger::new(y),
            z: SignedInteger::new(z),
            t: SignedInteger::new(t),
        },
        SUCCESS_RETVAL,
    )
}

pub fn iota_endomorphism<Fp2: Fp2Trait>(P: &Point<Fp2>) -> Point<Fp2> {
    let (xP, yP) = P.to_xy();

    Point::new_xy(&(-xP), &(Fp2::ZETA * yP))
}

pub fn frobenius_endomorphism<Fp2: Fp2Trait, const NUM_WORDS_P: usize>(
    p: &BigNum<NUM_WORDS_P>,
    P: &Point<Fp2>,
) -> Point<Fp2> {
    let (xP, yP) = P.to_xy();

    Point::new_xy(
        &xP.pow(&p.to_le_bytes(), p.nbits()),
        &yP.pow(&p.to_le_bytes(), p.nbits()),
    )
}

pub fn apply_endomorphism_from_quaternion<Fp2: Fp2Trait>(
    curve: &Curve<Fp2>,
    quat_basis: &[Point<Fp2>; 4],
    quat: &Quaternion,
) -> Point<Fp2> {
    let mut X = curve.mul(
        &quat_basis[0],
        &quat.x.n.to_digits::<u8>(Order::Lsf),
        quat.x.n.significant_bits() as usize,
    );
    X.set_condneg(quat.x.is_negative);
    let mut Y = curve.mul(
        &quat_basis[1],
        &quat.y.n.to_digits::<u8>(Order::Lsf),
        quat.y.n.significant_bits() as usize,
    );
    Y.set_condneg(quat.y.is_negative);
    let mut Z = curve.mul(
        &quat_basis[2],
        &quat.z.n.to_digits::<u8>(Order::Lsf),
        quat.z.n.significant_bits() as usize,
    );
    Z.set_condneg(quat.z.is_negative);
    let mut T = curve.mul(
        &quat_basis[3],
        &quat.t.n.to_digits::<u8>(Order::Lsf),
        quat.t.n.significant_bits() as usize,
    );
    T.set_condneg(quat.t.is_negative);

    let mut result = Point::INFINITY;
    curve.addto(&mut result, &X);
    curve.addto(&mut result, &Y);
    curve.addto(&mut result, &Z);
    curve.addto(&mut result, &T);

    result
}

/// Given a `basis` of the `p`^`e`-torsion and an `endo`morphism on `curve`, finds
/// the kernel of a degree-`p`^`e` isogeny backtracking the `endo`morphism
pub fn find_kernel_of_backtracking_isogeny_prime_power_degree<
    Fp2: Fp2Trait,
    const NUM_WORDS_DEG: usize,
>(
    curve: &Curve<Fp2>,
    quat_basis: &([Point<Fp2>; 4], [Point<Fp2>; 4]),
    endo: &Quaternion,
    reduced_degree: &BigNum<NUM_WORDS_DEG>,
) -> Point<Fp2> {
    // Apply the endomorphism to points basis points on the starting curve,
    // one of which will be used as the kernel for the dual endomorphism
    let endo_P = apply_endomorphism_from_quaternion(curve, &quat_basis.0, endo);
    let endo_Q = apply_endomorphism_from_quaternion(curve, &quat_basis.1, endo);

    // If endo(Q) is of full order, then it must be linearly independent to the kernel of
    // the backtracking isogeny we're trying to construct. Otherwise, endo(P) must be.
    let endoQ_has_full_order = !curve
        .mul(
            &endo_Q,
            &reduced_degree.to_le_bytes(),
            reduced_degree.nbits(),
        )
        .is_zero();
    let mut backtracking_isogeny_kernel = endo_P;
    backtracking_isogeny_kernel.set_cond(&endo_Q, endoQ_has_full_order);

    backtracking_isogeny_kernel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_mul_by_complex_power_test() {
        let test = |x: (&Integer, &Integer),
                    y: (&Integer, &Integer),
                    exp: u32,
                    cmp: (Integer, Integer)| {
            let result = complex_mul_by_complex_power(x, y, exp);
            assert_eq!(result, cmp);
        };

        test(
            (&Integer::from(3), &(-Integer::from(4))),
            (&Integer::from(1), &Integer::from(2)),
            0,
            (Integer::from(3), -Integer::from(4)),
        );
        test(
            (&Integer::from(3), &(-Integer::from(4))),
            (&Integer::from(1), &Integer::from(2)),
            1,
            (Integer::from(11), Integer::from(2)),
        );
        test(
            (&Integer::from(3), &(-Integer::from(4))),
            (&Integer::from(1), &Integer::from(2)),
            2,
            (Integer::from(7), Integer::from(24)),
        );
    }

    #[test]
    fn cornacchia_test() {
        let test = |n: Integer, p: Integer| {
            let (x, y, ok) = cornacchia(&n, &p);
            assert_eq!(ok, SUCCESS_RETVAL);
            assert!(x.square() + n.clone() * y.square() == p.clone());
        };

        test(Integer::from(1), Integer::from(2));
        test(Integer::from(1), Integer::from(5));
        test(Integer::from(1), Integer::from(41));
        test(Integer::from(3), Integer::from(7));
        test(Integer::from(1), Integer::from(29));
        test(Integer::from(1), Integer::from(1381));
    }

    #[test]
    fn cornacchia_test_fail() {
        let test = |n: Integer, p: Integer| {
            let (_, _, ok) = cornacchia(&n, &p);
            assert_eq!(ok, FAILURE_RETVAL);
        };

        test(Integer::from(3), Integer::from(2));
        test(Integer::from(5), Integer::from(5));
        test(Integer::from(82), Integer::from(41));
        test(Integer::from(21), Integer::from(7));
        test(Integer::from(29), Integer::from(29));
        test(Integer::from(2762), Integer::from(1381));
    }

    #[test]
    fn cornacchia_extended_prime_loop_test() {
        let test = |c: (&Integer, &Integer), p: &Integer, val: u32, cmp: (Integer, Integer)| {
            let (x, y, ok) = cornacchia_extended_prime_loop(c, p, val);
            assert_eq!(ok, SUCCESS_RETVAL);
            assert_eq!(x, cmp.0);
            assert_eq!(y, cmp.1);
        };

        test(
            (&Integer::from(1), &Integer::from(1)),
            &Integer::from(5),
            2,
            (-Integer::from(1), Integer::from(7)),
        );
        test(
            (&Integer::from(1), &Integer::from(0)),
            &Integer::from(2),
            2,
            (Integer::from(0), Integer::from(2)),
        );
        test(
            (&(-Integer::from(8)), &Integer::from(6)),
            &Integer::from(41),
            1,
            (-Integer::from(64), -Integer::from(2)),
        );
        test(
            (&Integer::from(31), &Integer::from(298)),
            &Integer::from(17),
            1,
            (-Integer::from(174), Integer::from(1223)),
        );

        let c_re = -Integer::from_str_radix("286292335164776146355279307454016", 10).unwrap();
        let c_im = Integer::from_str_radix("29390034047986289563880081110912", 10).unwrap();
        let cmp_re =
            Integer::from_str_radix("17365063775578586201598356324038842078208", 10).unwrap();
        let cmp_im =
            -Integer::from_str_radix("120145123962752358003585692097427531473856", 10).unwrap();
        test((&c_re, &c_im), &Integer::from(37), 11, (cmp_re, cmp_im));
    }

    #[test]
    fn cornacchia_extended_test() {
        let PRIMES = SMALL_PRIMES.map(|prime| Integer::from(prime));
        let BAD_PRIMES = SMALL_BAD_PRIMES.map(|prime| Integer::from(prime));

        let test = |n: Integer, primes: &[Integer]| {
            let (x, y, ok) = cornacchia_extended(&n, primes, &BAD_PRIMES);
            assert_eq!(ok, SUCCESS_RETVAL);
            assert_eq!(x.square() + y.square(), n);
        };

        let n = Integer::from_str_radix("158828126405", 10).unwrap();
        test(n, &PRIMES[..26]);

        let n = Integer::from_str_radix("40014512245953376820", 10).unwrap();
        test(n, &PRIMES[..26]);

        let n = Integer::from_str_radix("13588601138373462550928562093344019331858765165447160913867695800135179682916335386618687107041342641315567209049089762325285514045525398790223782031779675017956905089308664998960177516311434400605404622520285574320810358353904423423012810358448642090586663211232583668667454734852141666802903725775590093049043999776164077647933766930542685009426257611359117058369584270139317425716467676553177434794742587363282934876629123721404464944683384304420820986239594188922650094326930372850811146481429730669474269713207255284173968500838746261882707217336240321810667791730219009998376346822365457541240875528870506928821817510280272725782371744788162331496464288694956312256388500637865037597849249048688448700135956112736510275629520042877377629061438837268497413468138373466165224807808346472654909601360038557113093422495558353125000", 10).unwrap();
        test(n, &PRIMES[..100]);
    }
}
