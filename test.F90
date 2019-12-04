program foo

    use precise_sum

    implicit none


    real, dimension(5) :: density, q, n_re, nqd, qd

    real, parameter :: first_sum  = 1.48093748092651
    real, parameter :: second_sum = 1.48095517372973

    real :: correct
    real :: is, ks, ps, dcs


    density = (/1000.000, 1770.000, 1770.000, 1000.000, 1900.000/)
    n_re    = (/1.35934579372406, 1.48093748092651, 1.48093748092651, 1.52999997138977, 1.50999999046326/)

    correct = first_sum
    write (*,98) correct
    q       = (/0.0000000E+00, 1.1703869E-34, 1.0192583E-35, 0.0000000E+00,  0.0000000E+00/)

    !nqd = n_re * q/density
    !qd  =        q/density

    is  =       sum(n_re * q/density) /       sum(q/density)
    write (*,100) is, 100*abs(correct-is)/correct

    ks  = kahan_sum(n_re * q/density) / kahan_sum(q/density)
    write (*,101) ks, 100*abs(correct-ks)/correct

    ps  =  pair_sum(n_re * q/density) /  pair_sum(q/density)
    write (*,102) ps, 100*abs(correct-ps)/correct

    dcs = dcomp_sum(n_re * q/density) / dcomp_sum(q/density)
    write (*,103) dcs, 100*abs(correct-dcs)/correct


    write (*,*) ""
    correct = second_sum
    write (*,99) correct
    q       = (/0.0000000E+00, 6.0657997E-34, 2.0065308E-35, 0.0000000E+00,  4.0976052E-37/)

    is  =       sum(n_re * q/density) /       sum(q/density)
    write (*,100) is, 100*abs(correct-is)/correct

    ks  = kahan_sum(n_re * q/density) / kahan_sum(q/density)
    write (*,101) ks, 100*abs(correct-ks)/correct

    ps  =  pair_sum(n_re * q/density) /  pair_sum(q/density)
    write (*,102) ps, 100*abs(correct-ps)/correct

    dcs = dcomp_sum(n_re * q/density) / dcomp_sum(q/density)
    write (*,103) dcs, 100*abs(correct-dcs)/correct

 98 format (" First test, expected value: ", g30.12)
 99 format ("Second test, expected value: ", g30.12)
100 format ("Intrinsic sum: ",g30.12, " Percent diff: ", g30.12)
101 format ("    Kahan sum: ",g30.12, " Percent diff: ", g30.12)
102 format (" Pairwise sum: ",g30.12, " Percent diff: ", g30.12)
103 format ("Doub Comp sum: ",g30.12, " Percent diff: ", g30.12)

end program foo 
