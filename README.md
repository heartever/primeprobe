## prime+probe code targeting a given physical address

Basically our code for prime+probe cache attack is a simplified version of the code from the paper: Last-Level Cache Side-Channel Attacks are Practical.

Within the threat model of SGX the attacker can know the physical address of a target virtual address, as such in our code we provide a parameter for the targeted physical address (-a option). However you may need to write a kernel module to get the physical address for a target virtual address.

We have performed the attack against libgcrypt (the attached source files.) Also we verifed the attack on an unmodifed version of libgcrypt with the help of graphene-SGX. The target virtual address will be the address of the function mpih_sqr_n_basecase.
