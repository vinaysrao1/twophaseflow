# Brent's method (Brent 1973) for a bracketed scalar root.

"""
    brent(f, a, b; xtol=1e-12, rtol=4eps(), maxiter=200) -> x

Root of `f` in `[a, b]`; requires `f(a)` and `f(b)` of opposite sign.
"""
function brent(f, a::Float64, b::Float64; xtol=1e-12, rtol=4eps(), maxiter=200)
    fa, fb = f(a), f(b)
    fa == 0 && return a
    fb == 0 && return b
    sign(fa) == sign(fb) && throw(ArgumentError("root not bracketed: f($a)=$fa, f($b)=$fb"))
    c, fc = a, fa
    d = e = b - a
    for _ in 1:maxiter
        if sign(fb) == sign(fc)
            c, fc = a, fa
            d = e = b - a
        end
        if abs(fc) < abs(fb)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        end
        tol = 2rtol * abs(b) + xtol / 2
        m = (c - b) / 2
        (abs(m) <= tol || fb == 0) && return b
        if abs(e) >= tol && abs(fa) > abs(fb)
            s = fb / fa
            if a == c
                p, q = 2m * s, 1 - s
            else
                q, r = fa / fc, fb / fc
                p = s * (2m * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            end
            p > 0 ? (q = -q) : (p = -p)
            if 2p < min(3m * q - abs(tol * q), abs(e * q))
                e, d = d, p / q
            else
                d = e = m
            end
        else
            d = e = m
        end
        a, fa = b, fb
        b += abs(d) > tol ? d : copysign(tol, m)
        fb = f(b)
    end
    return b
end
