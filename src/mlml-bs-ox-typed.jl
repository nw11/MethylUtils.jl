using GSL

function expectation(M_B::Int16,u_B::Int16,p_m::Float64,p_h::Float64,betas::Vector{Float64})
	  log_beta = 0
    for k=0:M_B
        h_B = M_B - k
        m_B = k
    	  log_beta = sf_lnchoose(M_B, k) + m_B * log( p_m ) +  h_B * log ( p_h ) -  M_B * log( p_m + p_h )
        #println("k: $k, log_beta: $log_beta, p_m $p_m, p_h $p_h")
        betas[k+1] = exp(log_beta)
    end
    return betas
end

function maximisation(M_B::Int16, u_B::Int16, p_m::Float64, betas::Vector{Float64})
    (num,denom) = (0,0)
    #h_B = int16(0)
    #k=int16(0)
    for k=0:M_B
    	h_B = M_B - k
        #m_B = k
        num   +=  betas[k+1] * h_B * ( 1 - p_m)
        denom +=  betas[k+1] * ( u_B + h_B )
    end
    return num/denom
end

function update_pm(M_B::Int16, u_B::Int16, m_O::Int16, U_O::Int16, betas::Vector{Float64})
    num = m_O
    denom=0
    for k=0:M_B
    	m_B=k
    	num += betas[k+1]*m_B
    end
    denom = M_B + m_O + u_B + U_O
    return num/denom
end

 # Go for the estimate with the highest coverage??
function start_point(M_B,u_B,m_O,U_O)
        starting_estimate_p_m = M_B / (M_B + u_B)
        starting_estimate_p_h = 1 - starting_estimate_p_m
        return ( starting_estimate_p_m, starting_estimate_p_h )
end


# function returns the values
# p_m and p_h and the number of iterations taken to find the solution
function EM( M_B::Int16, u_B::Int16, m_O::Int16, U_O::Int16, p_m::Float64, p_h::Float64)
	tolerance = float64(1e-15)
	iter = int(0);
	verbose = 0
	betas = Array(Float64,M_B+1)
	while true
		betas = expectation(M_B, u_B, p_m, p_h,betas)
		p_m_old = float64(p_m)
    p_h_old = float64(p_h)
    if verbose == 1
       println("p_m: $p_m_old p_h: $p_h_old")
		end
		p_h   = maximisation(M_B, u_B, p_m, betas)
		p_m   = update_pm(M_B, u_B, m_O, U_O, betas)
		delta = max(abs(p_m_old - p_m), abs(p_h_old - p_h));
	      if iter > 500
            tolerance=delta
            if verbose == 1
                 println("Reached 500 iterations ( delta = $delta )")
            end
            return p_m,p_h, iter
        end

        if (delta <= tolerance)
        	if verbose == 1
                println( "delta <= tolerance ( $delta <= $tolerance ) - reached after $iter iterations ")
            end
            return p_m,p_h, iter
            #break
        end
        iter +=1
	end
end

function EM(BS,OX,count)

   if count % 10000 == 0
           t = TmStruct( time() )
           t_str = join( [t.hour,t.min,t.sec], ":")
           println("On $count. current time: ", t_str)
   end

   if any( x -> x < 0, (  BS.T_count, BS.C_count, OX.C_count, OX.T_count )  )
      println( "error negative counts: ", join( [BS.pos.chr,BS.pos.start, BS.T_count, BS.C_count, OX.C_count, OX.T_count], "\t" ) )
      return CytosineMethylation( BS.pos, BS.strand, 2, 2 )
   end

   if ! ( BS.pos.chr == OX.pos.chr && BS.pos.start == OX.pos.start )
     println("Error: BS and OX positions do not match up:  ", join( [BS.pos.chr,BS.pos.start, BS.T_count, BS.C_count, OX.pos.chr, OX.pos.start, OX.C_count, OX.T_count], "\t" ) )
   end

   est_p_m = ( OX.C_count / (OX.C_count + OX.T_count ))
   est_p_h = ( BS.C_count / ( BS.T_count + BS.C_count ) ) -  est_p_m

   if( est_p_h >= 0  )
        return CytosineMethylation( BS.pos, BS.strand, est_p_m, est_p_h )
   else
       (starting_estimate_p_m, starting_estimate_p_h) = start_point(BS.C_count,BS.T_count,OX.C_count,OX.T_count)
      (em_p_m,em_p_h) = EM(BS.C_count,BS.T_count,OX.C_count,OX.T_count,starting_estimate_p_m, starting_estimate_p_h)
       return CytosineMethylation( BS.pos, BS.strand, em_p_m, em_p_h )
    end
end

function run_mlml(BS,OX)
   println("Start run_mlml")
   return  pmap( EM, BS, OX, 1:length(BS) )
end

