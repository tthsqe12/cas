Ring interface: (no output x,y,z may alias any input a,b,c)

arithmetic:
    zero(x)         | x = 0
    one(x)          | x = 1
    set(x,a)        | x = a
    neg(x,a)        | x = -a
    neg(x,a)        | x = -x
    add(x,a,b)      | x = a + b
    add(x,a)        | x += a
    sub(x,a,b)      | x = a - b
    sub(x,a)        | x -= a
    mul(x,a,b)      | x = a*b
    mul(x,a)        | x *= a
    madd(x,a,b,c)   | x = a + b*c
    madd(x,a,b)     | x += a*b
    msub(x,a,b,c)   | x = a - b*c
    msub(x,a,b)     | x -= a*b
    pow(x,a,n)      | x = a^n
    inv(x,a)        | x = 1/a  with x*a == 1
    divexact(x,a,b) | x = a/b  with x*b == a
    divexact(x,a)   | x /= b   with x'*a == x

predicates:
