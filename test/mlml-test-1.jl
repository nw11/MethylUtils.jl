module mlmlTest

  using FactCheck
  using Distributions
  include(Pkg.dir("MethylUtils", "src","mlml-bs-ox-typed.jl"))
   include(Pkg.dir("MethylUtils", "src","Features.jl"))
  MLML_PATH = Pkg.dir("MethylUtils", "bin","mlml")


  function test_mlml(M_B, u_B, m_O, U_O, p_m,p_h)
      ratio_bs = M_B / (M_B+u_B)
      bs_total = M_B+u_B
      ratio_ox = m_O / (m_O+U_O)
      ox_total = m_O+U_O

      (em_p_m, em_p_h, iter) = EM(M_B, u_B, m_O, U_O, p_m,p_h)
      p_u = 1 - em_p_m - em_p_h
      println("EM: p_m: $em_p_m, p_h: $em_p_h, p_u: $p_u")

      lineA = ["chr1" "2" "3" "CpG:$bs_total" ratio_bs "+"]
      lineB = ["chr1" "2" "3" "CpG:$ox_total" ratio_ox "+"]
      println(typeof(lineA))
      writedlm( "bs.bed", lineA, "\t")
      writedlm( "ox.bed", lineB, "\t")

      cmd = `$MLML_PATH -tolerance 1e-15 -bsseq bs.bed -oxbsseq ox.bed`
      output = readall(cmd)
      (chr,start,stop,mc_txt,hmc,strand) = split(output, "\t")
      (text,mc) = split(mc_txt,":")
      ( mc,hmc) = map( parsefloat, [ mc,hmc])
      (em_p_m,em_p_h,mc,hmc) = map( x -> round(x,4),  [em_p_m,em_p_h,mc,hmc])
      return (em_p_m,em_p_h,mc,hmc)
  end


  function rand_CytosineCounts(num_cgs, max_cov, min_cov=0)
   coverage_distribution = [ x for x= 1:max_cov]
   cgs = Array(CytosineCount,num_cgs)
   k=1
   while k <= num_cgs
       (Cs, Ts) = rand_c_counts(0.8, coverage_distribution)
       if ( Cs + Ts ) >= min_cov
           cgs[k] = CytosineCount(GenomePos("chr1",1,2),true,int16(Cs) ,int16(Ts))
           k +=1
       end
   end
   return cgs
  end

  function rand_c_counts(prob, coverage_distribution)
    num_reads = sample(coverage_distribution)
    binom = Binomial(num_reads,prob)
    Cs = rand(binom)
    Ts = num_reads - Cs
    return (Cs,Ts)
  end

  # M_B = Counts of C's in Bis-Seq - including 5mc (m) and 5hmc (h)
  # u_B = counts of T's in Bis-Seq
  # m_O = Counts of C's in Oxidative Bis-Seq - 5mc
  # U_O = Counts of T's in Oxidative Bis-Seq - including 5hmc (h) and umc (u)

  (M_B,u_B,m_O,U_O) = map( int16, [5,5,7,3 ] )

  ##  M ratio  = 0.5 (5/10),  m ratio  = 0.7 (7/10) - overshoot example

  starting_estimate_p_m = M_B / (M_B + u_B)
  starting_estimate_p_h = 1 - starting_estimate_p_m
  betas = Array(Float64,M_B+1)
  betas = expectation(  M_B, u_B, starting_estimate_p_m, starting_estimate_p_h, betas )

  # beta has length M_B+1, as it is calculated over each possible combination of counts
  # for the latent variables (m_B, h_B) - i.e.  ( (0,5), (1,4), (2,3), ..., (5,0) )

  expected_betas = [0.03125,0.15624999999999994,0.3125,0.3125,0.15624999999999994,0.03125]

  # maximisation calculates p_h
  new_p_h = maximisation(M_B, u_B, starting_estimate_p_m, betas)
  new_p_m = update_pm(M_B, u_B, m_O, U_O, betas)
  (expected_new_p_h, expected_new_p_m) = (0.16666666666666666, 0.475)

  ( em_p_m, em_p_h, mlml_p_m, mlml_p_h) = test_mlml(M_B, u_B, m_O, U_O, starting_estimate_p_m,starting_estimate_p_h)

  println( "$em_p_m, $em_p_h, $mlml_p_m, $mlml_p_h")

  facts("expectation") do
   @fact length(betas) => M_B+1
   @fact betas => expected_betas
  end

  facts("maximisation") do
    @fact new_p_h => expected_new_p_h
    @fact new_p_m => expected_new_p_m
  end

  facts("EM1") do
    @fact em_p_m => mlml_p_m
    @fact em_p_h => mlml_p_h
  end

  # 10/15 (0.666)  10/20 (0.5) - No overshoot
  # M_B - m_O = p_h = 0.166
  (M_B,u_B,m_O,U_O) = map( int16, [10,5,10,10 ] )
  starting_estimate_p_m = M_B / (M_B + u_B)
  starting_estimate_p_h = 1 - starting_estimate_p_m

  #
  #
  expected_p_m = ( m_O / ( m_O + U_O ) )
  expected_p_h = starting_estimate_p_m - expected_p_m
  expected_p_h = round(expected_p_h,4)
  expected_p_m = round(expected_p_m,4)

  ( em_p_m, em_p_h, mlml_p_m, mlml_p_h) = test_mlml(M_B, u_B, m_O, U_O, starting_estimate_p_m,starting_estimate_p_h)
   println( "$em_p_m, $em_p_h, $mlml_p_m, $mlml_p_h")

 facts("EM2") do
    @fact em_p_m => expected_p_m
    @fact em_p_h => expected_p_h
  end


  # M_B = 90/190 (0.4736842)
  # m_O = 85/185 (0.4594595)
  # p_h = 0.01422475
  # No overshoot - but only by a small amount
  (M_B,u_B,m_O,U_O) = map( int16, [90,100,85,100 ] )
  starting_estimate_p_m = M_B / (M_B + u_B)
  starting_estimate_p_h = 1 - starting_estimate_p_m

  expected_p_m = ( m_O / ( m_O + U_O ) )
  expected_p_h = starting_estimate_p_m - expected_p_m
  expected_p_h = round(expected_p_h,4)
  expected_p_m = round(expected_p_m,4)

  ( em_p_m, em_p_h, mlml_p_m, mlml_p_h) = test_mlml(M_B, u_B, m_O, U_O, starting_estimate_p_m,starting_estimate_p_h)
   println( "$em_p_m, $em_p_h, $mlml_p_m, $mlml_p_h")

 facts("EM3") do
    @fact em_p_m => expected_p_m
    @fact em_p_h => expected_p_h
  end

  #### OTHER WAY AROUND ####
  # M_B = 85/185
  # m_O = 90/190
  # Overshoot - but only by a small amount
  (M_B,u_B,m_O,U_O) = map( int16, [85,100,90,100 ] )
  starting_estimate_p_m = M_B / (M_B + u_B)
  starting_estimate_p_h = 1 - starting_estimate_p_m

  expected_p_m = ( m_O / ( m_O + U_O ) )
  expected_p_h = starting_estimate_p_m - expected_p_m
  expected_p_h = round(expected_p_h,4)
  expected_p_m = round(expected_p_m,4)

  ( em_p_m, em_p_h, mlml_p_m, mlml_p_h) = test_mlml(M_B, u_B, m_O, U_O, starting_estimate_p_m,starting_estimate_p_h)
   println( "$em_p_m, $em_p_h, $mlml_p_m, $mlml_p_h")

  facts("EM3") do
    @fact em_p_m => mlml_p_m
    @fact em_p_h => mlml_p_h
  end

  # TEST MLML
  bs_counts = rand_CytosineCounts(10,100,5)
  ox_counts = rand_CytosineCounts(10,100,5)
  cgs = run_mlml( bs_counts , ox_counts )
  length(cgs)

  facts("mlml") do
    @fact length(cgs) => 10
  end

  # exitstatus()
end

# put in the other tests
# only calculate em when need to
# parallelise what I can
# put up to travis
