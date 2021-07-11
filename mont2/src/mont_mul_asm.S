.global fe_mont_mul_asm
.global _fe_mont_mul_asm

# See calculation in https://research.nccgroup.com/2021/06/09/optimizing-pairing-based-cryptography-montgomery-arithmetic-in-rust/
# See ../constant.py
.set NPRIME, 0x89f3fffcfffcfffd

.align 8
N:                                  # Field modulus for BLS12-381, LSB firsst
.quad  0xb9feffffffffaaab
.quad  0x1eabfffeb153ffff
.quad  0x6730d2a0f6b0f624
.quad  0x64774b84f38512bf
.quad  0x4b1ba7b6434bacd7
.quad  0x1a0111ea397fe69a

# Montgomery mult; assumes properly reduced input operands in 6x64-bit limbs
# %rdi holds address of result least significant limb
# %rsi and %rdx hold address of operands (a & b) least significant limb
_fe_mont_mul_asm:
fe_mont_mul_asm:

    lea     -55(%rip), %rcx         # Load address of field modulus into %rcx

    # Prologue: push all **modified** callee-save registers onto stack
    # We are free to use regs below and rax,rcx,rdx,rsi,rdi,r8,r9,r10,r11
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15
    pushq   %rbx

    # Note: %r10-%r15 hold limbs of working result; %r8 holds temp 'overflow'

    # First partial prod 'loop' derived from macro below, with incoming %r1x=0
    xorq    %r8, %r8                # Clear flags and save zero for use later
    movq    %rdx, %r9               # Hang on to address of b.v[0]
    movq    (%r9), %rdx             # Load b.v[0] into %rdx for following mulx
    mulxq   0(%rsi), %r10, %rbx     # a.v[0] * %rdx:b.v[0] -> lo:%r10, hi:%rbx
    mulxq   8(%rsi), %r11, %rax     # a.v[1] * %rdx:b.v[0] -> lo:%r11, hi:%rax
    adcxq   %rbx, %r11              # Add earlier hi into later lo as %r11
    mulxq   16(%rsi), %r12, %rbx    # a.v[2] * %rdx:b.v[0] -> lo:%r12, hi:%rbx
    adcxq   %rax, %r12              # Add earlier hi into later lo as %r12
    mulxq   24(%rsi), %r13, %rax    # a.v[3] * %rdx:b.v[0] -> lo:%r13, hi:%rax
    adcxq   %rbx, %r13              # Add earlier hi into later lo as %r13
    mulxq   32(%rsi), %r14, %rbx    # a.v[4] * %rdx:b.v[0] -> lo:%r13, hi:%rbx
    adcxq   %rax, %r14              # Add earlier hi into later lo as %r14
    mulxq   40(%rsi), %r15, %rax    # a.v[5] * %rdx:b.v[0] -> lo:%r15, hi:%r8
    adcxq   %rbx, %r15              # Add earlier hi into later lo as %r15
    adcxq   %rax, %r8               # Propogate prior carry into %rbp

    # calculate m and drop it into %rdx for use in subsequent mulx
    movq    $NPRIME, %rdx
    imul    %r10, %rdx              # Puts least significant 64b into %rdx

    # First reduction step, base address of N[0-6] is %rcx
    xorq    %rax, %rax              # Clear flags
    mulxq   0(%rcx), %rax, %rbx     # N[0] * %rdx:m -> lo:%rax, %hi:rbx
    adcxq   %r10, %rax              # %rax discarded, but generate carry out
    adoxq   %rbx, %r11              # partial_a[0]
    mulxq   8(%rcx), %r10, %rbx     # N[1] * %rdx:m -> lo:%r10, %hi:rbx
    adcxq   %r11, %r10              # A[0] in %r10 for next round
    adoxq   %rbx, %r12              # partial_a[1]
    mulxq   16(%rcx), %r11, %rbx    # N[2] * %rdx:m -> lo:%r11, %hi:rbx
    adcxq   %r12, %r11              # A[1] in %r11 for next round
    adoxq   %rbx, %r13              # partial_a[2]
    mulxq   24(%rcx), %r12, %rbx    # N[3] * %rdx:m -> lo:%r12, %hi:rbx
    adcxq   %r13, %r12              # A[2] in %r12 for next round
    adoxq   %rbx, %r14              # partial_a[3]
    mulxq   32(%rcx), %r13, %rbx    # N[4] * %rdx:m -> lo:%r13, %hi:rbx
    adcxq   %r14, %r13              # A[3] in %r13 for next round
    adoxq   %rbx, %r15              # partial_a[4]
    mulxq   40(%rcx), %r14, %rbx    # N[5] * %rdx:m -> lo:%r14, %hi:rbx
    adcxq   %r15, %r14              # A[4] in %r14 for next round
    movq    $0, %r15                # Clear %r15; need to sum two carry_in
    adcxq   %r8, %r15               # partial_a[5]
    adoxq   %rbx, %r15              # A[5] in %r15 for next round

# The 'inner loop' as a macro, to be instantiated 5 times w/ specific offset
.macro partial_product offset:req
    xorq    %r8, %r8                # Clear flags and save zero for use later
    movq    \offset(%r9), %rdx      # Load b.v[i] into %rdx for following mulx
    mulxq   0(%rsi), %rax, %rbx     # a.v[0] * %rdx:b.v[i] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r10              # Add lo into %r10 for red[0]
    adoxq   %rbx, %r11              # Add hi into %r11 as partial_red[1]
    mulxq   8(%rsi), %rax, %rbx     # a.v[1] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r11              # Add lo into %r11 as red[1]
    adoxq   %rbx, %r12              # Add hi into %r12 as partial_red[2]
    mulxq   16(%rsi), %rax, %rbx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r12              # Add lo into %r11 as red[1]
    adoxq   %rbx, %r13              # Add hi into %r13 as partial_red[3]
    mulxq   24(%rsi), %rax, %rbx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r13              # Add lo into %r13 as red[2]
    adoxq   %rbx, %r14              # Add hi into %r14 as partial_red[4]
    mulxq   32(%rsi), %rax, %rbx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r14              # Add lo into %r14 as red[4]
    adoxq   %rbx, %r15              # Add hi into %r15 as partial_red[5]
    mulxq   40(%rsi), %rax, %rbx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r15              # Add lo into %r15 as red[5]
    adoxq   %r8, %rbx               # Add hi into %rbp as partial_red[6]
    adcxq   %rbx, %r8               # bring carry_in to red[6]

    # calculate m and drop it into %rdx
    movq    $NPRIME, %rdx
    imul    %r10, %rdx

    # Reduction step, base address of N[0-6] is %rcx
    xorq    %rax, %rax              # Clear flags
    mulxq   0(%rcx), %rax, %rbx     # N[0] * %rdx:m -> lo:%rax, %hi:rbx
    adcxq   %r10, %rax              # %rax discarded, but generate carry out
    adoxq   %rbx, %r11              # partial_a[0]
    mulxq   8(%rcx), %r10, %rbx     # N[1] * %rdx:m -> lo:%r10, %hi:rbx
    adcxq   %r11, %r10              # A[0] in %r10 for next round
    adoxq   %rbx, %r12              # partial_a[1]
    mulxq   16(%rcx), %r11, %rbx    # N[2] * %rdx:m -> lo:%r11, %hi:rbx
    adcxq   %r12, %r11              # A[1] in %r11 for next round
    adoxq   %rbx, %r13              # partial_a[2]
    mulxq   24(%rcx), %r12, %rbx    # N[3] * %rdx:m -> lo:%r12, %hi:rbx
    adcxq   %r13, %r12              # A[2] in %r12 for next round
    adoxq   %rbx, %r14              # partial_a[3]
    mulxq   32(%rcx), %r13, %rbx    # N[4] * %rdx:m -> lo:%r13, %hi:rbx
    adcxq   %r14, %r13              # A[3] in %r13 for next round
    adoxq   %rbx, %r15              # partial_a[4]
    mulxq   40(%rcx), %r14, %rbx    # N[5] * %rdx:m -> lo:%r14, %hi:rbx
    adcxq   %r15, %r14              # A[4] in %r14 for next round
    movq    $0, %r15                # Clear %r15; need to sum two carry_in
    adcxq   %r8, %r15               # partial_a[5]
    adoxq   %rbx, %r15              # A[5] in %r15 for next round
.endm

    # Instantiate 'inner loop' 5 times (unroll b.v[1-5])
    partial_product 8
    partial_product 16
    partial_product 24
    partial_product 32
    partial_product 40

    # Make a copy of the result to prepare for subtracting modulus
    movq    %r10, %r8
    movq    %r11, %r9
    movq    %r12, %rax
    movq    %r13, %rbx
    movq    %r14, %rdx
    movq    %r15, %rsi

    # Subtract the modulus
    subq    0(%rcx), %r8
    sbbq    8(%rcx), %r9
    sbbq    16(%rcx), %rax
    sbbq    24(%rcx), %rbx
    sbbq    32(%rcx), %rdx
    sbbq    40(%rcx), %rsi

    # If there was final borrow, we must store original
    cmovcq  %r10, %r8               # Unlikely to be constant-time...
    cmovcq  %r11, %r9               # ...but program flow is 'same'
    cmovcq  %r12, %rax
    cmovcq  %r13, %rbx
    cmovcq  %r14, %rdx
    cmovcq  %r15, %rsi

    # Store result
    movq    %r8, 0(%rdi)
    movq    %r9, 8(%rdi)
    movq    %rax, 16(%rdi)
    movq    %rbx, 24(%rdi)
    movq    %rdx, 32(%rdi)
    movq    %rsi, 40(%rdi)

    # Epilogue: pop all callee-save registers from the stack
    popq    %rbx
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    ret
