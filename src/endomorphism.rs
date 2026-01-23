//! Functions for transforming quaternions to endomorphisms and associated helper functions.
//!
//! Adapted from the KLaPoTi-Rust library: https://github.com/isogeny-klapoti/klapoti-rust.

use fp2::traits::Fp2 as Fp2Trait;
use isogeny::elliptic::{basis::BasisX, curve::Curve, projective_point::Point};
use rug::{
    Integer,
    integer::{IsPrime, Order},
    rand::RandState,
};

use crate::{FAILURE_RETVAL, SUCCESS_RETVAL, bn::BigNum};

// First 1006 primes
const SMALL_PRIMES: [u32; 1006] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
    937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039,
    1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279,
    1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409,
    1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
    1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
    1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
    1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
    1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113,
    2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251,
    2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477,
    2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647,
    2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731,
    2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857,
    2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001,
    3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163,
    3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299,
    3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407,
    3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539,
    3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659,
    3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919,
    3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051,
    4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201,
    4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
    4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463,
    4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603,
    4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903,
    4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009,
    5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153,
    5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303,
    5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441,
    5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
    5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843,
    5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987,
    6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131,
    6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269,
    6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373,
    6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553,
    6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691,
    6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829,
    6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967,
    6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109,
    7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247,
    7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451,
    7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559,
    7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
    7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841,
    7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963,
];
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

    let target_norm = Integer::from_digits(&target_norm.to_le_bytes(), Order::Lsf);
    let p = Integer::from_digits(&p.to_le_bytes(), Order::Lsf);

    // FIXME: does this do rounding division? And does rounding division mess up our parameters?
    let m_prime = (&target_norm / p.clone()).sqrt();
    let m_prime_range = &TWO * m_prime.clone() + Integer::ONE;

    let mut x = Integer::ZERO;
    let mut y = Integer::ZERO;
    let mut z = Integer::ZERO;
    let mut t = Integer::ZERO;

    let mut counter = 1_000_000;
    let mut found_solution = FAILURE_RETVAL;
    while counter > 0 && found_solution == FAILURE_RETVAL {
        z = m_prime_range.clone().random_below(&mut rng) - &m_prime;

        let m_dblprime = (&target_norm / p.clone() - z.clone().square()).sqrt();
        let m_dblprime_range = &TWO * m_dblprime.clone() + Integer::ONE;
        t = m_dblprime_range.clone().random_below(&mut rng) - &m_dblprime;

        let norm = &target_norm - &p * (z.clone().square() + t.clone().square());
        if norm.is_probably_prime(30) == IsPrime::No {
            counter -= 1;
            continue;
        }

        (x, y, found_solution) = cornacchia(Integer::ONE, &norm);
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

pub fn apply_endomorphism_from_quaternion<Fp2: Fp2Trait, const NUM_WORDS_P: usize>(
    p: &BigNum<NUM_WORDS_P>,
    curve: &Curve<Fp2>,
    quat: &Quaternion,
    P: &Point<Fp2>,
) -> Point<Fp2> {
    // TODO: precompute this for 2^a basis points
    fn i_endo<Fp2: Fp2Trait>(P: &Point<Fp2>) -> Point<Fp2> {
        let (xP, yP) = P.to_xy();

        Point::new_xy(&(-xP), &(Fp2::ZETA * yP))
    }

    // TODO: precompute this for 2^a basis points
    fn frob_endo<Fp2: Fp2Trait, const NUM_WORDS_P: usize>(
        p: &BigNum<NUM_WORDS_P>,
        P: &Point<Fp2>,
    ) -> Point<Fp2> {
        let (xP, yP) = P.to_xy();

        Point::new_xy(
            &xP.pow(&p.to_le_bytes(), p.nbits()),
            &yP.pow(&p.to_le_bytes(), p.nbits()),
        )
    }

    let mut X = curve.mul(
        P,
        &quat.x.n.to_digits::<u8>(Order::Lsf),
        quat.x.n.significant_bits() as usize,
    );
    X.set_condneg(quat.x.is_negative);
    let mut Y = curve.mul(
        &i_endo(P),
        &quat.y.n.to_digits::<u8>(Order::Lsf),
        quat.y.n.significant_bits() as usize,
    );
    Y.set_condneg(quat.y.is_negative);
    let mut Z = curve.mul(
        &frob_endo(p, P),
        &quat.z.n.to_digits::<u8>(Order::Lsf),
        quat.z.n.significant_bits() as usize,
    );
    Z.set_condneg(quat.z.is_negative);
    let mut T = curve.mul(
        &i_endo(&frob_endo(p, P)),
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
    const NUM_WORDS_P: usize,
    const NUM_WORDS_DEG: usize,
>(
    field_characteristic: &BigNum<NUM_WORDS_P>,
    curve: &Curve<Fp2>,
    endo: &Quaternion,
    basis: &BasisX<Fp2>,
    reduced_degree: &BigNum<NUM_WORDS_DEG>,
) -> Point<Fp2> {
    // Apply the endomorphism to points basis points on the starting curve,
    // one of which will be used as the kernel for the dual endomorphism
    let (P, Q) = curve.lift_basis(basis);
    let endo_P = apply_endomorphism_from_quaternion(field_characteristic, curve, endo, &P);
    let endo_Q = apply_endomorphism_from_quaternion(field_characteristic, curve, endo, &Q);

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
