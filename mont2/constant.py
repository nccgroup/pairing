# python3 constant.py

N = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 2**384

def extended_euclidean(a, b):
    previous_x, x = 1, 0
    previous_y, y = 0, 1
    while b:
        q = a//b
        x, previous_x = previous_x - q*x, x
        y, previous_y = previous_y - q*y, y
        a, b = b, a % b
    return a, previous_x, previous_y

a, Rp, Np = extended_euclidean(R, N)         # R*Rp + N*Np == 1
Rp1 = Rp + N; Np1 = -1 * (Np - R)            # Rp1 and Np1 are positive
print(a == 1, 0 < Rp1 < N, 0 < Np1 < R, R*Rp1 - N*Np1 == 1)
print(f'Rp1 = {hex(Rp1)}\nNp1 = {hex(Np1)}')
print(N*Np1 == R*Rp1 - 1, N*Np1 % R == R-1)  # N*Np1 = -1 mod R
