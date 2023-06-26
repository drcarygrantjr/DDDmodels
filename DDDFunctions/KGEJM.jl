"""
    kge(q_sim, q_obs)
Compute modified Kling-Gupta efficiency
by Jan Magnusson
"""

function kge(q_sim, q_obs)

    # ikeep = q_obs .!= -999.

    if all(isnan, q_sim) || all(isnan, q_obs)

        kge = NaN

    else

        ikeep = .!isnan.(q_obs)

        q_sim .= q_sim[ikeep]
        q_obs .= q_obs[ikeep]

        r = cor(q_sim, q_obs)

        beta = mean(q_sim) / mean(q_obs)

        gamma = (std(q_sim) / mean(q_sim)) / (std(q_obs) / mean(q_obs))

        kge = 1 - sqrt( (r-1)^2 + (beta-1)^2 + (gamma-1)^2 )

    end

    return kge

end


