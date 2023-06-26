"""
    nse(q_sim, q_obs)
Compute Nash-Sutcliffe efficiency
by Jan Magnusson
"""

function nse(q_sim, q_obs)

    # ikeep = q_obs .!= -999.

    ikeep = .!isnan.(q_obs)

    q_sim .= q_sim[ikeep]
    q_obs .= q_obs[ikeep]

    ns = 1 - sum((q_sim .- q_obs).^2) / sum((q_obs .- mean(q_obs)).^2)

end

