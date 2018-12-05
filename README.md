## prime+probe code targeting a given physical address

Basically our code for prime+probe cache attack is a simplified version of the code from the paper: Last-Level Cache Side-Channel Attacks are Practical.

Within the threat model of SGX the attacker can know the physical address of a target virtual address, as such in our code we provide a parameter for the targeted physical address (-a option). However you may need to write a kernel module to get the physical address for a target virtual address.

We have performed the attack against libgcrypt (the attached source files.) Also we verifed the attack on an unmodifed version of libgcrypt with the help of graphene-SGX. The target virtual address will be the address of the function mpih_sqr_n_basecase.

Attention should be paid to the cache slice mapping which is hardcoded in function [getslicemapping(...)](https://github.com/heartever/primeprobe/blob/0178e80793ddb18f00e2c106a6a8e91d548e6dd4/newattack.c#L201) of newattack.c. Our testbed is equipped with an i7-6700k processor with hyperthreading enabled. It is likely to work for other processors with 4 physical cores when hyperthreading is turned on.

If hyperthreading is turned off or not supported, you may try other configurations provided in the paper (e.g. the first pair of h1 and h0). I am not sure about the slice mapping for processors with 6 cores supporting SGX.
