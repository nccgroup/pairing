.global fe_mont_mul_asm
.global _fe_mont_mul_asm

# Modulus from least significant limb towards most significant
.set N0, 0xb9feffffffffaaab
.set N1, 0x1eabfffeb153ffff
.set N2, 0x6730d2a0f6b0f624
.set N3, 0x64774b84f38512bf
.set N4, 0x4b1ba7b6434bacd7
.set N5, 0x1a0111ea397fe69a

# See calculation in https://research.nccgroup.com/2021/06/09/optimizing-pairing-based-cryptography-montgomery-arithmetic-in-rust/
.set NPRIME, 0x89f3fffcfffcfffd

# Montgomery multiplication; assumes properly reduced operands in six 64-bit limbs
# %rdi holds address of result least significant limb
# %rsi and %rdx hold address of each operand (a & b) least significant limb
_fe_mont_mul_asm:
fe_mont_mul_asm:

    # Prologue: push all **modified** callee-save registers onto stack
    pushq   %r12
    pushq   %r13
    pushq   %r14
    pushq   %r15

    # Note: %r10-%r15 hold the limbs of our working result; %r8 holds temp 'overflow'

    # First partial prod code derived from macro below, assumes incoming %r1x=0
    xorq    %r8, %r8                # Clear flags and save zero for use later
    movq    %rdx, %r9               # Hang on to address of b.v[0]
    movq    (%r9), %rdx             # Load value b.v[0] into %rdx for mulx
    mulxq   0(%rsi), %r10, %rcx     # a.v[0] * %rdx:b.v[0] -> lo:%r10, hi:%rcx
    mulxq   8(%rsi), %r11, %rax     # a.v[1] * %rdx:b.v[0] -> lo:%r11, hi:%rax
    adcxq    %rcx, %r11              # Add earlier hi into later lo as %r11
    mulxq   16(%rsi), %r12, %rcx    # a.v[2] * %rdx:b.v[0] -> lo:%r12, hi:%rcx
    adcxq    %rax, %r12              # Add earlier hi into later lo as %r12
    mulxq   24(%rsi), %r13, %rax    # a.v[3] * %rdx:b.v[0] -> lo:%r13, hi:%rax
    adcxq    %rcx, %r13              # Add earlier hi into later lo as %r13
    mulxq   32(%rsi), %r14, %rcx    # a.v[4] * %rdx:b.v[0] -> lo:%r13, hi:%rcx
    adcxq    %rax, %r14              # Add earlier hi into later lo as %r14
    mulxq   40(%rsi), %r15, %rax     # a.v[5] * %rdx:b.v[0] -> lo:%r15, hi:%r8
    adcxq    %rcx, %r15              # Add earlier hi into later lo as %r15
    adcxq    %rax, %r8                 # Propogate prior carry into %rbp

    # calculate m and drop it into %rdx for use in subsequent mulx
    movq    $NPRIME, %rdx
    imul    %r10, %rdx              # Saves least significant 64b into %rdx

    # First reduction step, base address of MOD.v[0-6] was passed in as %rcx
    xorq    %rax, %rax                # Clear flags and save zero for use later
    movq    $N0, %rcx
    mulxq   %rcx, %rax, %rcx     # MOD.v[0] * %rdx:m -> lo:%rax, %hi:rbx
    adcxq   %r10, %rax              # %rax discarded, but generate carry out
    adoxq   %rcx, %r11              # partial_a[0]
    movq    $N1, %rcx
    mulxq   %rcx, %r10, %rcx     # MOD.v[1] * %rdx:m -> lo:%r10, %hi:rbx
    adcxq   %r11, %r10              # A[0] in %r10 for next round
    adoxq   %rcx, %r12              # partial_a[1]
    movq    $N2, %rcx
    mulxq   %rcx, %r11, %rcx    # MOD.v[2] * %rdx:m -> lo:%r11, %hi:rbx
    adcxq   %r12, %r11              # A[1] in %r11 for next round
    adoxq   %rcx, %r13              # partial_a[2]
    movq    $N3, %rcx
    mulxq   %rcx, %r12, %rcx    # MOD.v[3] * %rdx:m -> lo:%r12, %hi:rbx
    adcxq   %r13, %r12              # A[2] in %r12 for next round
    adoxq   %rcx, %r14              # partial_a[3]
    movq    $N4, %rcx
    mulxq   %rcx, %r13, %rcx    # MOD.v[4] * %rdx:m -> lo:%r13, %hi:rbx
    adcxq   %r14, %r13              # A[3] in %r13 for next round
    adoxq   %rcx, %r15              # partial_a[4]
    movq    $N5, %rcx
    mulxq   %rcx, %r14, %rcx    # MOD.v[5] * %rdx:m -> lo:%r14, %hi:rbx
    adcxq   %r15, %r14              # A[4] in %r14 for next round
    movq    $0, %r15               # r8; clear %r15 because we need two carry_in
    adcxq   %r8, %r15              # partial_a[5]
    adoxq   %rcx, %r15              # A[5] in %r15 for next round



# The inner loop as a macro, to be instantiated 6 times
# %rsi holds address of a.[v0], **%rbx** holds address of b.v[i]
.macro partial_product offset:req
    xorq    %r8, %r8
    movq    \offset(%r9), %rdx      # Load b.v[i] into %rdx for implicit mulx
    mulxq   0(%rsi), %rax, %rcx     # a.v[0] * %rdx:b.v[i] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r10              # Add lo into %r10 for red[0]
    adoxq   %rcx, %r11              # Add hi into %r11 as partial_red[1]
    mulxq   8(%rsi), %rax, %rcx     # a.v[1] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r11              # Add lo into %r11 as red[1]
    adoxq   %rcx, %r12              # Add hi into %r12 as partial_red[2]
    mulxq   16(%rsi), %rax, %rcx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r12              # Add lo into %r11 as red[1]
    adoxq   %rcx, %r13              # Add hi into %r13 as partial_red[3]
    mulxq   24(%rsi), %rax, %rcx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r13              # Add lo into %r13 as red[2]
    adoxq   %rcx, %r14              # Add hi into %r14 as partial_red[4]
    mulxq   32(%rsi), %rax, %rcx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r14              # Add lo into %r14 as red[4]
    adoxq   %rcx, %r15              # Add hi into %r15 as partial_red[5]
    mulxq   40(%rsi), %rax, %rcx    # a.v[2] * %rdx:b.v[1] -> lo:%rax, hi:%rbx
    adcxq   %rax, %r15              # Add lo into %r15 as red[5]
    adoxq   %r8, %rcx               # Add hi into %rbp as partial_red[6]
    adcxq   %rcx, %r8               # bring carry_in to red[6]

    # calculate m and drop it into %rdx
    movq    $NPRIME, %rdx
    imul    %r10, %rdx

    # Put address of modulus into %rsi
    xorq    %rax, %rax
    movq    $N0, %rcx
    mulxq   %rcx, %rax, %rcx
    adcxq   %r10, %rax
    adoxq   %rcx, %r11
    movq    $N1, %rcx
    mulxq   %rcx, %r10, %rcx
    adcxq   %r11, %r10
    adoxq   %rcx, %r12
    movq    $N2, %rcx
    mulxq   %rcx, %r11, %rcx
    adcxq   %r12, %r11
    adoxq   %rcx, %r13
    movq    $N3, %rcx
    mulxq   %rcx, %r12, %rcx
    adcxq   %r13, %r12
    adoxq   %rcx, %r14
    movq    $N4, %rcx
    mulxq   %rcx, %r13, %rcx
    adcxq   %r14, %r13
    adoxq   %rcx, %r15
    movq    $N5, %rcx
    mulxq   %rcx, %r14, %rcx
    adcxq   %r15, %r14
    movq    $0, %r15
    adoxq   %rcx, %r15
    adcxq   %r8, %r15
.endm

    partial_product 8
    partial_product 16
    partial_product 24
    partial_product 32
    partial_product 40

    # Make a copy of the result
    movq  %r10, %r8
    movq  %r11, %r9
    movq  %r12, %rax
    movq  %r13, %rcx
    movq  %r14, %rdx
    movq  %r15, %rsi

    # Subtract the modulus; we are short one register, so push/pop rdi
    pushq   %rdi
    movq    $N0, %rdi
    subq    %rdi, %r8        # subtract LSB of modulus  ### possibly recalculate r8 again to save a register?
    movq    $N1, %rdi
    sbbq    %rdi, %r9
    movq    $N2, %rdi
    sbbq    %rdi, %rax
    movq    $N3, %rdi
    sbbq    %rdi, %rcx
    movq    $N4, %rdi
    sbbq    %rdi, %rdx
    movq    $N5, %rdi
    sbbq    %rdi, %rsi      # subtract MSB of modulus  TODO: <--- this is where we run out of registers!!
    popq    %rdi

    # If there was no carry above, then prepare to store that
    cmovcq  %r10, %r8
    cmovcq  %r11, %r9
    cmovcq  %r12, %rax
    cmovcq  %r13, %rcx
    cmovcq  %r14, %rdx
    cmovcq  %r15, %rsi

    # Store result
    movq    %r8, 0(%rdi)
    movq    %r9, 8(%rdi)
    movq    %rax, 16(%rdi)
    movq    %rcx, 24(%rdi)
    movq    %rdx, 32(%rdi)
    movq    %rsi, 40(%rdi)

    # Epilogue: pop all callee-save registers from the stack
    popq    %r15
    popq    %r14
    popq    %r13
    popq    %r12
    #popq    %rbx
    ret
