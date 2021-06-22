## Modeling the natural history of human diseases
#' Compute the iPSE in natural history of human diseases
#' @param data a data.frame where columns must be in the order of T1, d1, T2, d2, T3, d3, Z
#' @return A list of self-described result
#' @export
DNHiPSE = function(data, zvec = NULL, plot_result = TRUE, iPSEweights = TRUE, weights = NULL, bootstrap = FALSE, boot_times = 1000, num_of_cores = 1, for_bootstrap = FALSE,
                   push_warnings = TRUE, plot_unit = 1, save_plot = FALSE, plot_name = NULL, folder_name = NULL, match_ylim = FALSE, PM = FALSE, timer = TRUE, Med_TTEM = FALSE, sensitivity_analysis = FALSE, sensitivity_analysis_match_ylim = FALSE){
  # zvec = NULL; plot_result = TRUE; iPSEweights = TRUE; weights = NULL; bootstrap = TRUE; boot_times = 1000; num_of_cores = 1; for_bootstrap = FALSE; push_warnings = TRUE; plot_unit = 1; save_plot = F; plot_name = NULL; folder_name = NULL; match_ylim = FALSE; PM = FALSE; timer = TRUE; Med_TTEM = TRUE; sensitivity_analysis = TRUE; sensitivity_analysis_match_ylim = 'effect'
  if(bootstrap & num_of_cores > 1){require(snow); require(doSNOW); require(pracma); require(foreach);}

  ## data pre-check
  if(!for_bootstrap){
    pre_check = DNH_pre_check_data(data, push_warnings = push_warnings)
    data = pre_check$data
    org_data = pre_check$org_data
    min_diff = pre_check$min_diff
  }
  if(is.vector(data)){return(list(index_wrong_order = data))}

  ## get components
  if(is.null(weights)){weights = rep(1, dim(data)[1])}
  if(!is.vector(weights)){stop("weights must be a numeric vector.")}
  components = DNH_get_components(data, weights)

  ## need plug-in
  zvec_len = ifelse(is.null(zvec), 4, length(zvec))
  effect = DNH_get_effect(components, zvec, iPSEweights)
  if(!(PM == FALSE)){proportion_mediated = DNH_PM(effect, PM)}

  # only give estimators when T3 is observed
  if(!for_bootstrap){for(i in 1:zvec_len){effect[[i]] = effect[[i]][components$observed, ]}}

  # get MED_TTEM
  if(Med_TTEM & (!for_bootstrap)){MED_effect = DNH_Med_TTEM(components)}

  # get sensitivity analysis
  if(sensitivity_analysis & (!for_bootstrap)){
    sen_ana = DNH_get_sen_ana(data, components, zvec, iPSEweights, names(effect))
    for(i in 1:length(sen_ana)){
      for(j in 1:length(sen_ana[[1]])){
        sen_ana[[i]][[j]] = sen_ana[[i]][[j]][components$observed, ]
      }
    }
  }
  ## bootstrap
  if(bootstrap){
    if(!(PM == FALSE)){PM_list = list(rep(0, boot_times), rep(0, boot_times), rep(0, boot_times), rep(0, boot_times))}

    if(num_of_cores > 1){
      boot_effect_mat_iPSE = vector(mode = 'list', length = zvec_len)
      for(i in 1:zvec_len){boot_effect_mat_iPSE[[i]] = matrix(0, nrow = boot_times, ncol = sum(data$d3))}
      cl = snow::makeCluster(num_of_cores[1])
      my_functions = c('DNH_simulation', 'DNH_generate_data', 'DNH_get_wb_true', 'DNH_get_wa_true', 'DNH_simulation_true_value',
                       'DNH_change_Z_to_12', 'DNH_boot_shift_data', "DNH_pre_check_data",
                       'DNH_leave', 'DNH_make_plot_main', 'DNH_list2names', 'DNH_rep_row', 'DNH_revcumsum', 'DNH_sort_mat', 'DNH_nan2zero', 'DNH_get_position', 'DNH_get_sick_number',
                       'DNH_get_Ybar', 'DNH_get_Ybar_n1..', 'DNH_get_Ybar_..n2', 'DNH_get_Ybar_n1n2', 'DNH_get_w_n1n2', 'DNH_get_w_n1..', 'DNH_get_w_..n2', 'DNH_get_wcond_n1n2',
                       'DNH_get_weights', 'DNH_get_dN3bar_n1n2', 'DNH_get_wa', 'DNH_get_wb', 'DNH_get_dL_n1n2', 'DNH_get_components', 'DNH_get_counterfactual_hazard',
                       'DNH_get_effect', 'DNH_PM', 'DNHiPSE')
      snow::clusterExport(cl, my_functions); doSNOW::registerDoSNOW(cl); pb = txtProgressBar(max = boot_times, style = 3); progress = function(n) setTxtProgressBar(pb, n); opts = list(progress = progress)

      boot_effect = foreach(i = 1:boot_times, .options.snow = opts, .combine = 'c') %dopar%{
        set.seed(20210303 + i)
        boot_index = sample(1:dim(data)[1], dim(data)[1], replace = TRUE)
        boot_data = data[boot_index, ]
        boot_data = DNH_boot_shift_data(boot_data, min_diff)
        boot_data = boot_data[sort(boot_data$T3, index.return = TRUE)$ix, ]
        boot_effect = tryCatch(DNHiPSE(boot_data, zvec = zvec, iPSEweights = iPSEweights, weights = weights, PM = PM, bootstrap = FALSE, plot_result = FALSE, for_bootstrap = TRUE), error = function(msg){return(list(msg = msg))})
        rm(list = "boot_data")
        gc()
        return(list(boot_effect))
      }

      snow::stopCluster(cl); pracma::fprintf('\n');
      index_safe_boot = NULL
      for(i in 1:boot_times){
        if(is.null(boot_effect[[i]]$msg)){
          for(j in 1:zvec_len){
            boot_effect_mat_iPSE[[j]][i, ] = approx(boot_effect[[i]]$effect[[j]]$time, boot_effect[[i]]$effect[[j]]$effect, effect[[1]]$time, method = 'linear', yleft = 0, rule = 2, ties = 'max')$y
            if(!(PM == FALSE)){PM_list[[j]][i] = boot_effect[[i]]$PM[j, 1]}
          }
          index_safe_boot = c(index_safe_boot, i)
        }
      }
      for(i in 1:zvec_len){
        boot_effect_mat_iPSE[[i]] = DNH_sort_mat(boot_effect_mat_iPSE[[i]][index_safe_boot, ])
        effect[[i]]$boot_lower = boot_effect_mat_iPSE[[i]][round(dim(boot_effect_mat_iPSE[[i]])[1] * 0.025), ]
        effect[[i]]$boot_upper = boot_effect_mat_iPSE[[i]][round(dim(boot_effect_mat_iPSE[[i]])[1] * 0.975), ]
      }

    }else{
      safe_num = 0
      inner_counter = 0
      index_safe_boot = NULL
      boot_effect_mat_iPSE = vector(mode = 'list', length = zvec_len)
      for(i in 1:zvec_len){boot_effect_mat_iPSE[[i]] = matrix(0, nrow = 25 * 2, ncol = sum(data$d3))}
      boot_effect_mat_iPSE_tmp = boot_effect_mat_iPSE

      boot_times = (round(boot_times/50)) * 50
      boot_times = ifelse(boot_times <= 50, 50, boot_times)
      if(timer){pracma::fprintf('| bootstrap        20        30        40        50        60        70        80        90        | 100\n')}
      for(i in 1:boot_times){
        if(timer & (i %% (boot_times/100) == 0)){pracma::fprintf('-')}
        set.seed(20210303 + i)
        boot_index = sample(1:dim(data)[1], dim(data)[1], replace = TRUE)
        boot_data = data[boot_index, ]
        boot_data = DNH_boot_shift_data(boot_data, min_diff)
        boot_data = boot_data[sort(boot_data$T3, index.return = TRUE)$ix, ]
        boot_effect = tryCatch(DNHiPSE(boot_data, zvec = zvec, iPSEweights = iPSEweights, weights = weights, PM = PM, bootstrap = FALSE, plot_result = FALSE, for_bootstrap = TRUE), error = function(msg){return(list(msg = msg))})
        if(is.null(boot_effect$msg)){
          index_safe_boot = c(index_safe_boot, i)
          safe_num = safe_num + 1
          if(safe_num <= 50){
            for(j in 1:zvec_len){
              if(!(PM == FALSE)){PM_list[[j]][i] = boot_effect$PM[j, 1]}
              boot_effect_mat_iPSE[[j]][safe_num, ] = approx(boot_effect$effect[[j]]$time, boot_effect$effect[[j]]$effect, effect[[1]]$time, method = 'linear', yleft = 0, rule = 2, ties = 'max')$y
            }
          }else{
            inner_counter = inner_counter + 1
            for(j in 1:zvec_len){
              if(!(PM == FALSE)){PM_list[[j]][i] = boot_effect$PM[j, 1]}
              boot_effect_mat_iPSE_tmp[[j]][inner_counter, ] = approx(boot_effect$effect[[j]]$time, boot_effect$effect[[j]]$effect, effect[[1]]$time, method = 'linear', yleft = 0, rule = 2, ties = 'max')$y
            }
            if(inner_counter == 50){
              inner_counter = 0
              for(j in 1:zvec_len){
                boot_mat_tmp = rbind(boot_effect_mat_iPSE[[j]], boot_effect_mat_iPSE_tmp[[j]])
                boot_mat_tmp = DNH_sort_mat(boot_mat_tmp)
                boot_effect_mat_iPSE[[j]] = boot_mat_tmp[c(1:25, 76:100), ]
              }
            }
          }
        }else{
          warning(boot_effect$msg)
        }
      }
      for(i in 1:zvec_len){
        effect[[i]]$boot_lower = boot_effect_mat_iPSE[[i]][25, ]
        effect[[i]]$boot_upper = boot_effect_mat_iPSE[[i]][26, ]
      }
    }

    if(!(PM == FALSE)){
      upper_PM = data.frame(rep(0, zvec_len), row.names = row.names(proportion_mediated))
      lower_PM = data.frame(rep(0, zvec_len), row.names = row.names(proportion_mediated))
      for(i in 1:zvec_len){
        PM_list[[i]] = sort(PM_list[[i]])
        upper_PM[i, 1] = PM_list[[i]][boot_times * 0.975]
        lower_PM[i, 1] = PM_list[[i]][boot_times * 0.025]
      }
      proportion_mediated = data.frame(PM = proportion_mediated, upper_PM = upper_PM, lower_PM = lower_PM)
      colnames(proportion_mediated) = c('PM', 'upper_PM', 'lower_PM')
    }

    if(push_warnings & (length(index_safe_boot) != boot_times)){
      index_wrong_boot = paste(setdiff(1:boot_times, index_safe_boot), collapse = ", ")
      warning("Bootstrap number", index_wrong_boot, "have some bugs. Please tell me.")
    }
  }

  if(plot_result & bootstrap){
    if(!save_plot){DNH_plot_result(effect, zvec, plot_unit, match_ylim)}
    if(save_plot){DNH_save_plot_result(effect, zvec, plot_unit, plot_name, folder_name, match_ylim)}
  }
  if(sensitivity_analysis & plot_result & (!for_bootstrap)){
    if(!save_plot){DNH_plot_sen_ana(sen_ana, zvec, plot_unit, sensitivity_analysis_match_ylim)}
    if(save_plot){DNH_save_plot_sen_ana(sen_ana, zvec, plot_unit, plot_name, folder_name, sensitivity_analysis_match_ylim)}
  }

  if(!for_bootstrap){
    result = list(effect = effect, data = org_data, components = components)
  }else{
    result = list(effect = effect)
  }

  if(!(PM == FALSE)){result$PM = proportion_mediated}
  if(Med_TTEM & (!for_bootstrap)){result$Med_TTEM = MED_effect}
  if(sensitivity_analysis & (!for_bootstrap)){result$sensitivity_analysis = sen_ana}
  tmp = DNH_leave()
  return(result)
}

# generate sample data and simulation
#' @export
DNH_simulation = function(type, hypo, sample_size = 1e3, weights = NULL, iPSEweights = TRUE, repeat_size = 1e3, num_of_cores = 1, ylim = NULL){
  if(hypo == 'null'){
    alphaZ = 0; betaZ = 0; gammaZ = 0
  }else{
    alphaZ = 1; betaZ = 1; gammaZ = 1
  }

  if(type == "coverage"){
    timeline = seq(0, 10, by = 0.01)
    true_value = DNH_simulation_true_value(alphaZ, betaZ, gammaZ, out_time = timeline, iPSEweights = iPSEweights)

    if(repeat_size <= 200){warning("Small repeat_size may lead to ill result while computing coverage", immediate. = TRUE)}
    coverage_table = list(rep(0, length(timeline)), rep(0, length(timeline)), rep(0, length(timeline)), rep(0, length(timeline)))

    if(num_of_cores > 1){
      cl = snow::makeCluster(num_of_cores[1])
      my_functions = c('DNH_simulation', 'DNH_generate_data', 'DNH_get_wb_true', 'DNH_get_wa_true', 'DNH_simulation_true_value',
                       'DNH_change_Z_to_12', 'DNH_boot_shift_data', "DNH_pre_check_data",
                       'DNH_leave', 'DNH_make_plot_main', 'DNH_list2names', 'DNH_rep_row', 'DNH_revcumsum', 'DNH_sort_mat', 'DNH_nan2zero', 'DNH_get_position', 'DNH_get_sick_number',
                       'DNH_get_Ybar', 'DNH_get_Ybar_n1..', 'DNH_get_Ybar_..n2', 'DNH_get_Ybar_n1n2', 'DNH_get_w_n1n2', 'DNH_get_w_n1..', 'DNH_get_w_..n2', 'DNH_get_wcond_n1n2',
                       'DNH_get_weights', 'DNH_get_dN3bar_n1n2', 'DNH_get_wa', 'DNH_get_wb', 'DNH_get_dL_n1n2', 'DNH_get_components', 'DNH_get_counterfactual_hazard',
                       'DNH_get_effect', 'DNH_PM', 'DNHiPSE')
      snow::clusterExport(cl, my_functions); doSNOW::registerDoSNOW(cl); pb = txtProgressBar(max = repeat_size, style = 3); progress = function(n) setTxtProgressBar(pb, n); opts = list(progress = progress)

      result_tmp = foreach(i = 1:repeat_size, .options.snow = opts, .combine = 'c') %dopar%{
        data = DNH_generate_data(m = sample_size, alphaZ, betaZ, gammaZ, seed = i)
        result_tmp = tryCatch(DNHiPSE(data, bootstrap = TRUE, iPSEweights = iPSEweights, weights = weights, plot_result = FALSE)$effect, error = function(msg){return(list(msg = msg))})

        boot_lower_now = vector("list", 4)
        boot_upper_now = vector("list", 4)
        if(is.null(result_tmp$msg)){
          for(j in 1:4){
            boot_lower_now[[j]] = approx(result_tmp[[j]]$time, result_tmp[[j]]$boot_lower, timeline, yleft = 0, rule = 2)$y
            boot_upper_now[[j]] = approx(result_tmp[[j]]$time, result_tmp[[j]]$boot_upper, timeline, yleft = 0, rule = 2)$y
          }
        }
        rm(list = 'data')
        gc()
        return(list(boot_lower_now, boot_upper_now))
      }
      snow::stopCluster(cl); pracma::fprintf('\n');
      for(j in 1:4){
        for(i in 1:repeat_size){coverage_table[[j]] = coverage_table[[j]] + ((result_tmp[[i * 2 - 1]][[j]] - true_value[[j]]$effect) * (result_tmp[[2 * i]][[j]] - true_value[[j]]$effect) <= 0)}
        coverage_table[[j]] = coverage_table[[j]]/repeat_size
      }
    }else{
      for(i in 1:repeat_size){
        data = DNH_generate_data(m = sample_size, alphaZ, betaZ, gammaZ, seed = i)
        result_tmp = DNHiPSE(data, bootstrap = TRUE, iPSEweights = iPSEweights, weights = weights, plot_result = FALSE, timer = FALSE)$effect

        for(j in 1:4){
          boot_lower_now = approx(result_tmp[[j]]$time, result_tmp[[j]]$boot_lower, timeline, yleft = 0, rule = 2)$y
          boot_upper_now = approx(result_tmp[[j]]$time, result_tmp[[j]]$boot_upper, timeline, yleft = 0, rule = 2)$y
          coverage_table[[j]] = coverage_table[[j]] + ((boot_lower_now - true_value[[j]]$effect) * (boot_upper_now - true_value[[j]]$effect) <= 0)
        }
      }
    }
    data_true = DNH_generate_data(m = 1e6, alphaZ, betaZ, gammaZ, seed = 1)
    observed_index = data_true$d1 & data_true$d2 & data_true$d3
    look_at_index = c(sort(data_true$T1[observed_index])[sum(observed_index) * 0.75], 0, sort(data_true$T3[observed_index])[sum(observed_index) * 0.25]) * 100
    look_at_index[2] = (look_at_index[1] + look_at_index[3])/2

    coverage = list(percentile = data.frame(matrix(0, nrow = 4, ncol = 3), row.names = c("ZY", "Z2Y", "Z1Y", "Z12Y")))
    colnames(coverage$percentile) = c('follow_up_25', 'follow_up_50', 'follow_up_75')
    for(j in 1:4){coverage$percentile[j, ] = coverage_table[[j]][look_at_index]}
    # coverage$raw = data.frame(time = timeline, coverage = coverage_table)
    # colnames(coverage$raw) = c("time", "ZY", "Z2Y", "Z1Y", "Z12Y")
    return(coverage)

  }else{
    data_true = DNH_generate_data(m = 1e6, alphaZ, betaZ, gammaZ, seed = 1)
    observed_index = data_true$d1 & data_true$d2 & data_true$d3
    max_time = sort(data_true$T3[observed_index])[sum(observed_index) * 0.75] * 1.05
    timeline = seq(0, max_time, by = 0.01)
    true_value = DNH_simulation_true_value(alphaZ, betaZ, gammaZ, out_time = timeline, iPSEweights = iPSEweights)
    if(is.null(ylim)){
      ylimvec = list(c(0, 0), c(0, 0), c(0, 0), c(0, 0))
      ylim_all = c(0, 0)
    }
    pracma::fprintf('| unbiasedness     20        30        40        50        60        70        80        90        | 100\n')
    result = lapply(seq_len(repeat_size), function(x){matrix(0, nrow = length(timeline), ncol = 4)})
    for(i in 1:repeat_size){
      if(i %% (repeat_size/100) == 0){pracma::fprintf("-")}
      data = DNH_generate_data(m = sample_size, alphaZ, betaZ, gammaZ, seed = i)
      result_tmp = DNHiPSE(data, bootstrap = FALSE, iPSEweights = iPSEweights, weights = weights)$effect

      for(j in 1:4){
        result[[i]][, j] = approx(result_tmp[[j]]$time, result_tmp[[j]]$effect, timeline, yleft = 0, rule = 2)$y
        if(is.null(ylim)){
          ylimvec[[j]][1] = min(ylimvec[[j]][1], result[[i]][, j][result[[i]][, j] < max_time])
          ylimvec[[j]][2] = max(ylimvec[[j]][2], result[[i]][, j][result[[i]][, j] < max_time])
          ylim_all[1] = min(ylim_all[1], ylimvec[[j]][1])
          ylim_all[2] = max(ylim_all[2], ylimvec[[j]][2])
        }
      }
    }

    mean_vec = list(rep(0, length(timeline)), rep(0, length(timeline)), rep(0, length(timeline)), rep(0, length(timeline)))
    for(j in 1:4){
      main = DNH_make_plot_main(j)
      if(is.null(ylim)){plot(NULL, xlim = c(0, max_time), ylim = ylimvec[[j]], xlab = 'time', ylab = 'effect', main = main)}
      if(!is.null(ylim)){plot(NULL, xlim = c(0, max_time), ylim = ylim[[j]], xlab = 'time', ylab = 'effect', main = main)}

      # plot(NULL, xlim = c(0, max_time), ylim = ylim_all, xlab = 'time', ylab = 'effect', main = main)
      for(i in 1:repeat_size){
        mean_vec[[j]] = mean_vec[[j]] + result[[i]][, j]
        lines(timeline, result[[i]][, j], lwd = 0.3)
      }
      mean_vec[[j]] = mean_vec[[j]]/repeat_size
      lines(true_value[[j]], col = 'green', lwd = 2)
      lines(timeline, mean_vec[[j]], col = 'red', lwd = 2)
      legend("bottomleft", c("average", "true value"), fill = c(2, 3), bty = 'n')
      abline(h = 0, col = 'darkgrey')
    }
    if(is.null(ylim)){ylim = ylimvec}
    return(list(ylim = ylim))
  }
}
#' @export
DNH_generate_data = function(m, alphaZ, betaZ, gammaZ, seed = 1996, censoring = TRUE, underlying = FALSE, messy_data = FALSE, T1toT2 = 1, T1toT3 = 1, T2toT3 = 1){
  set.seed(seed+25)
  Z = c(rep(0, m/2), rep(1, m/2))

  trueT1 = rweibull(m, shape = 1, scale = exp(alphaZ * Z))
  trueT2 = 0.5 * rweibull(m, shape = 1, scale = exp(betaZ * Z)) + T1toT2 * trueT1
  trueT3 = 0.5 * rweibull(m, shape = 1, scale = exp(gammaZ * Z)) + T1toT3 * trueT1 + T2toT3 * trueT2

  if(messy_data){
    trueT1 = floor(trueT1 * 100)/100
    trueT2 = floor(trueT2 * 100)/100
    trueT3 = floor(trueT3 * 100)/100
  }
  C = pmin(rweibull(m, shape = 2, scale = 10), runif(m, 10, 15)) + 100 * (!censoring)

  T3 = pmin(trueT3, C)
  T2 = pmin(trueT2, trueT3, C)
  T1 = pmin(trueT1, trueT2, trueT3, C)

  d3 = (T3 == trueT3)
  d2 = (T2 == trueT2)
  d1 = (T1 == trueT1)

  if(underlying){
    observedData = data.frame(T1, d1, T2, d2, T3, d3, Z)
    underlyingData = data.frame(trueT1, trueT2, trueT3, C, Z)
    return(list(observedData = observedData, underlyingData = underlyingData))
  }else{
    return(data.frame(T1, d1, T2, d2, T3, d3, Z))
  }
}

# true value
#' @export
DNH_get_wb_true = function(wcondn1n2, wn1.., zvec){
  zb = zvec[2]; zc = zvec[3]; zd = zvec[4];
  if(FALSE){
    wb00 = wcondn1n2$wcond00[[zb]] * wn1..$w0.[[zc]]
    wb01 = wcondn1n2$wcond10[[zb]] * wn1..$w0.[[zc]]
    wb10 = wcondn1n2$wcond01[[zb]] * wn1..$w1.[[zc]]
    wb11 = wcondn1n2$wcond11[[zb]] * wn1..$w1.[[zc]]
  }else{
    wb00 = (wcondn1n2$wcond00[[zb]] * wn1..$w0.[[zd]] + wcondn1n2$wcond01[[zb]] * wn1..$w1.[[zd]]) * wn1..$w0.[[zc]]
    wb01 = (wcondn1n2$wcond10[[zb]] * wn1..$w0.[[zd]] + wcondn1n2$wcond11[[zb]] * wn1..$w1.[[zd]]) * wn1..$w0.[[zc]]
    wb10 = (wcondn1n2$wcond00[[zb]] * wn1..$w0.[[zd]] + wcondn1n2$wcond01[[zb]] * wn1..$w1.[[zd]]) * wn1..$w1.[[zc]]
    wb11 = (wcondn1n2$wcond10[[zb]] * wn1..$w1.[[zd]] + wcondn1n2$wcond11[[zb]] * wn1..$w1.[[zd]]) * wn1..$w1.[[zc]]
  }

  wbn1n2 = list(wb00 = wb00, wb01 = wb01, wb10 = wb10, wb11 = wb11)
  return(wbn1n2)
}
#' @export
DNH_get_wa_true = function(weights, wbn1n2, zvec){
  za = zvec[1]
  wa00 = 1/(wbn1n2$wb00 + wbn1n2$wb01 + wbn1n2$wb10 + wbn1n2$wb11 * weights[[za]])
  wa01 = wa00
  wa10 = wa00
  wa11 = weights[[za]]/(wbn1n2$wb00 + wbn1n2$wb01 + wbn1n2$wb10 + wbn1n2$wb11 * weights[[za]])

  wan1n2 = list(wa00 = wa00, wa01 = wa01, wa10 = wa10, wa11 = wa11)
  return(wan1n2)
}
#' @export
DNH_simulation_true_value = function(alphaZ, betaZ, gammaZ, iPSEweights = TRUE, out_time = NULL, onlydL11 = FALSE){
  # out_time = time_pt
  ## Z is zero and one
  z0 = 0
  z1 = 1

  if(is.null(out_time)){
    time_pt = seq(0, 10, by = 1e-3)
  }else{
    time_pt = out_time
  }

  parType = 1
  g1 = list(Z0 = NULL, Z1 = NULL); g2 = g1; g3 = g1;
  g1$Z0 = (parType == 1) * exp(-alphaZ * z0) + (parType == 2) * (1 + alphaZ * z0)
  g1$Z1 = (parType == 1) * exp(-alphaZ * z1) + (parType == 2) * (1 + alphaZ * z1)

  g2$Z0 = (parType == 1) * exp(-betaZ * z0) * 2 + (parType == 2) * (1 + betaZ * z0)
  g2$Z1 = (parType == 1) * exp(-betaZ * z1) * 2 + (parType == 2) * (1 + betaZ * z1)

  g3$Z0 = (parType == 1) * exp(-gammaZ * z0) * 2 + (parType == 2) * (1 + gammaZ * z0)
  g3$Z1 = (parType == 1) * exp(-gammaZ * z1) * 2 + (parType == 2) * (1 + gammaZ * z1)

  ## compute survival function
  survival1 = list(Z0 = NULL, Z1 = NULL); survival2 = survival1; survival3 = survival1; density3 = list(Z0 = NULL, Z1 = NULL);
  survival1$Z0 = 1 - sdprisk::phypoexp(time_pt, g1$Z0)
  survival1$Z1 = 1 - sdprisk::phypoexp(time_pt, g1$Z1)

  survival2$Z0 = 1 - sdprisk::phypoexp(time_pt, c(g1$Z0, g2$Z0 * (1 + 1e-6)))
  survival2$Z1 = 1 - sdprisk::phypoexp(time_pt, c(g1$Z1, g2$Z1 * (1 + 1e-6)))

  density3$Z0 = sdprisk::dhypoexp(time_pt, c(g1$Z0/2, g2$Z0 * (1 + 1e-6), g3$Z0 * (1 + 2e-6)))
  density3$Z1 = sdprisk::dhypoexp(time_pt, c(g1$Z1/2, g2$Z1 * (1 + 1e-6), g3$Z1 * (1 + 2e-6)))
  survival3$Z0 = 1 - sdprisk::phypoexp(time_pt, c(g1$Z0/2, g2$Z0 * (1 + 1e-6), g3$Z0 * (1 + 2e-6)))
  survival3$Z1 = 1 - sdprisk::phypoexp(time_pt, c(g1$Z1/2, g2$Z1 * (1 + 1e-6), g3$Z1 * (1 + 2e-6)))

  ## true wn1..
  w0. = list(Z0 = survival1$Z0/survival3$Z0, Z1 = survival1$Z1/survival3$Z1)
  w1. = list(Z0 = 1 - w0.$Z0, Z1 = 1 - w0.$Z1)
  # plot(time_pt, w0.$Z0, type = 'l', col = 'red'); lines(components$wn1..$w0.$Z1)
  # plot(time_pt, w1.$Z0, type = 'l', col = 'red'); lines(components$wn1..$w1.$Z1)
  wn1.. = list(w0. = w0., w1. = w1.)

  ## true wcond
  wcond00 = list(Z0 = 1, Z1 = 1)
  wcond01 = list(Z0 = (survival2$Z0 - survival1$Z0)/(survival3$Z0 - survival1$Z0), Z1 = (survival2$Z1 - survival1$Z1)/(survival3$Z1 - survival1$Z1))
  wcond01$Z0[1] = 1; wcond01$Z1[1] = 1
  wcond10 = list(Z0 = 0, Z1 = 0)
  wcond11 = list(Z0 = 1 - wcond01$Z0, Z1 = 1 - wcond01$Z1)
  # plot(time_pt, rep(1, length(time_pt)), type = 'l', col = 'red'); lines(components$wcondn1n2$wcond00$Z1)
  # plot(time_pt, rep(0, length(time_pt)), type = 'l', col = 'red'); lines(components$wcondn1n2$wcond10$Z1)
  # plot(time_pt, wcond01$Z0, type = 'l', col = 'red'); lines(components$wcondn1n2$wcond01$Z1)
  # plot(time_pt, wcond11$Z0, type = 'l', col = 'red'); lines(components$wcondn1n2$wcond11$Z1)
  wcondn1n2 = list(wcond00 = wcond00, wcond01 = wcond01, wcond10 = wcond10, wcond11 = wcond11)

  ## true dL
  dL11 = list(Z0 = NULL, Z1 = NULL)
  dL11[[1]] = density3[[1]]/(survival3[[1]] - survival2[[1]]); dL11[[1]][1] = 0.5; dL11[[1]] = dL11[[1]] * c(0, diff(time_pt))
  dL11[[2]] = density3[[2]]/(survival3[[2]] - survival2[[2]]); dL11[[2]][1] = 0.5; dL11[[2]] = dL11[[2]] * c(0, diff(time_pt))
  weights = list(Z0 = NULL, Z1 = NULL)
  if(iPSEweights){
    weights[[1]] = exp(-cumsum(dL11[[1]]))
    weights[[2]] = exp(-cumsum(dL11[[2]]))
  }else{
    weights[[1]] = rep(1, length(time_pt))
    weights[[2]] = rep(1, length(time_pt))
  }

  # plot(time_pt, weights$Z0, type = 'l', col = 'red'); lines(components$iPSEweights$weights11$Z1)
  if(onlydL11){return(data.frame(time_pt, dL11))}

  ## compute ipse
  zveclist = list(c(2, 1, 1, 1), c(1, 1, 1, 1),
                  c(2, 2, 1, 1), c(2, 1, 1, 1),
                  c(2, 2, 2, 1), c(2, 2, 1, 1),
                  c(2, 2, 2, 2), c(2, 2, 2, 1))
  effect_tmp = vector(mode = 'list', length = 8)
  for(i in 1:8){
    zvec = zveclist[[i]] ## za etc are served as index. Therefore, they are 0/1 or 1/2 doesn't matter.
    # plot(time_pt, DNH_get_wb_true(wcondn1n2, wn1.., zvec)$wb00, type = 'l', col = 'red')
    # lines(wbn1n2$wb00$time, wbn1n2$wb00$value)
    # plot(time_pt, DNH_get_wb_true(wcondn1n2, wn1.., zvec)$wb01, type = 'l', col = 'red')
    # lines(wbn1n2$wb01$time, wbn1n2$wb01$value)
    # plot(time_pt, DNH_get_wb_true(wcondn1n2, wn1.., zvec)$wb10, type = 'l', col = 'red')
    # lines(wbn1n2$wb10$time, wbn1n2$wb10$value)
    # plot(time_pt, DNH_get_wb_true(wcondn1n2, wn1.., zvec)$wb11, type = 'l', col = 'red')
    # lines(wbn1n2$wb11$time, wbn1n2$wb11$value)
    wbn1n2 = DNH_get_wb_true(wcondn1n2, wn1.., zvec)
    wan1n2 = DNH_get_wa_true(weights, wbn1n2, zvec)

    za = zvec[1]
    effect_tmp[[i]] = cumsum(wan1n2$wa11 * wbn1n2$wb11 * dL11[[za]])
  }
  # plot(time_pt, effect_tmp[[i]], type = 'l', col = 'red')
  # lines(effect_tmp1[[i]]$time, effect_tmp1[[i]]$value)

  effect = list(PSE1000 = NULL, PSE1100 = NULL, PSE1110 = NULL, PSE1111 = NULL)
  effect$PSE1000 = data.frame(time = time_pt, effect = effect_tmp[[1]] - effect_tmp[[2]])
  effect$PSE1100 = data.frame(time = time_pt, effect = effect_tmp[[3]] - effect_tmp[[4]])
  effect$PSE1110 = data.frame(time = time_pt, effect = effect_tmp[[5]] - effect_tmp[[6]])
  effect$PSE1111 = data.frame(time = time_pt, effect = effect_tmp[[7]] - effect_tmp[[8]])

  if(!is.null(out_time)){
    effect$PSE1000 = data.frame(time = out_time, effect = approx(effect$PSE1000$time, effect$PSE1000$effect, out_time, method = 'linear', rule = 2, yleft = 0)$y)
    effect$PSE1100 = data.frame(time = out_time, effect = approx(effect$PSE1100$time, effect$PSE1100$effect, out_time, method = 'linear', rule = 2, yleft = 0)$y)
    effect$PSE1110 = data.frame(time = out_time, effect = approx(effect$PSE1110$time, effect$PSE1110$effect, out_time, method = 'linear', rule = 2, yleft = 0)$y)
    effect$PSE1111 = data.frame(time = out_time, effect = approx(effect$PSE1111$time, effect$PSE1111$effect, out_time, method = 'linear', rule = 2, yleft = 0)$y)
  }
  return(effect)
}

# data preprocess
#' @export
DNH_pre_check_data = function(data, push_warnings = TRUE){
  data = as.data.frame(data)
  colnames(data) = c('T1', 'd1', 'T2', 'd2', 'T3', 'd3', 'Z')
  ## check data type
  if(sum(is.na(data))){
    NA_columns = which(apply(data, 2, function(x){sum(is.na(x))}) > 0)
    stop("NA exist(s) in column: ", paste(NA_columns, collapse = ', '), ".")
  }
  if(sum(data$d1 %in% c(1, 0)) != dim(data)[1]){stop("data$d1 should be logical or 0/1.")}
  if(sum(data$d2 %in% c(1, 0)) != dim(data)[1]){stop("data$d2 should be logical or 0/1.")}
  if(sum(data$d3 %in% c(1, 0)) != dim(data)[1]){stop("data$d3 should be logical or 0/1.")}
  if(dim(data)[2] > 7){stop("The current method cannot adjust confounder. Please remove them.")}
  if(dim(data)[2] < 7){stop("The number of columns of data should be 7 and in the order of T1, d1, T2, d2, T3, d3, and Z.")}
  if(sum(data$d1) == 0){stop("No T1 are observed.")}
  if(sum(data$d2) == 0){stop("No T2 are observed.")}
  if(sum(data$d3) == 0){stop("No T3 are observed.")}
  unique_Z = sort(unique(data$Z))
  if(length(unique_Z) != 2){stop("data$Z should contain exactly two distinct elements.")}
  data$Z = DNH_change_Z_to_12(data$Z, unique_Z)

  data$d1 = (data$d1 == 1)
  data$d2 = (data$d2 == 1)
  data$d3 = (data$d3 == 1)

  ## check weird index
  ## negative data
  negative_index = sort(unique(c(which(data$T1 < 0), which(data$T2 < 0), which(data$T3 < 0))))
  if(length(negative_index) > 0){
    if(push_warnings){warning("The following number of rows have negative time: ", paste(negative_index, collapse = ', '), ".\nThey will be removed.")}
    data = data[-negative_index, ]
  }

  ## check logic
  T1_less_than_T2_but_unobserved = which((data$T1 < data$T2) & (data$d1 == 0))
  T1_less_than_T3_but_unobserved = which((data$T1 < data$T3) & (data$d1 == 0) & (data$d2 == 0))
  T2_less_than_T3_but_unobserved = which((data$T2 < data$T3) & (data$d2 == 0))
  if(sum(T1_less_than_T2_but_unobserved) + sum(T1_less_than_T3_but_unobserved) + sum(T2_less_than_T3_but_unobserved) > 0){
    index_wrong_order = paste(sort(unique(c(T1_less_than_T2_but_unobserved, T1_less_than_T3_but_unobserved, T2_less_than_T3_but_unobserved))), collapse = ', ')
    if(push_warnings){warning("When T1 is unobserved, T1 must = T2 (or T3). When T2 is unobserved, T2 must = T3.\nThe following number of rows should be checked again: ", index_wrong_order, ".")}
    data$d1[T1_less_than_T2_but_unobserved] = TRUE
    data$d1[T1_less_than_T3_but_unobserved] = TRUE
    data$d2[T2_less_than_T3_but_unobserved] = TRUE
  }

  ## get minimal difference
  min_diff = min(diff(unique(sort(c(data$T1, data$T2, data$T3))))) / 1.5

  ## > and < relation
  orderT1T2 = which(data$d1 == 1)[data$T1[data$d1 == 1] > data$T2[data$d1 == 1]]
  orderT1T3 = which(data$d1 == 1)[data$T1[data$d1 == 1] > data$T3[data$d1 == 1]]
  orderT2T3 = which(data$d2 == 1)[data$T2[data$d2 == 1] > data$T3[data$d2 == 1]]
  if(sum(orderT1T2) + sum(orderT1T3) + sum(orderT2T3) > 0){
    index_wrong_order = paste(sort(unique(c(orderT1T2, orderT1T3, orderT2T3))), collapse = ', ')
    if(push_warnings){warning("When T1 is observed, T1 must <= T2 (or T3). When T2 is observed, T2 must <= T3.\nThe following number of rows should be checked again: ", index_wrong_order, ".")}
    return(index_wrong_order)
  }

  ## check ties
  diffT1 = diff(sort(data$T1)); tieT1 = sum(diffT1 == 0);
  diffT2 = diff(sort(data$T2)); tieT2 = sum(diffT2 == 0);
  diffT3 = diff(sort(data$T3)); tieT3 = sum(diffT3 == 0);
  if(push_warnings & tieT1 > 0){warning("Ties are detected in T1. I will move them by a small amount. Use $data to fetch the actual data.")}
  if(push_warnings & tieT2 > 0){warning("Ties are detected in T2. I will move them by a small amount. Use $data to fetch the actual data.")}
  if(push_warnings & tieT3 > 0){warning("Ties are detected in T3. I will move them by a small amount. Use $data to fetch the actual data.")}
  if(tieT1 + tieT2 + tieT3 > 0){
    min_diff = min(diffT1[diffT1 != 0], diffT2[diffT2 != 0], diffT3[diffT3 != 0])/2
    shift = (1:dim(data)[1]) / dim(data)[1] * min_diff
    shift = sample(shift, length(shift))
    data$T1 = data$T1 + shift; data$T2 = data$T2 + shift; data$T3 = data$T3 + shift;
  }

  ## check if T1 and T2 are the same (when observed)
  sameT1T2 = which(data$T1 == data$T2 & data$d1 == 1 & data$d2 == 1)
  sameT1T3 = which(data$T1 == data$T3 & data$d1 == 1 & data$d3 == 1)
  sameT2T3 = which(data$T2 == data$T3 & data$d2 == 1 & data$d3 == 1)
  if(sum(sameT1T2) + sum(sameT1T3) + sum(sameT2T3) > 0){
    shift = min(diffT1[diffT1 != 0], diffT2[diffT2 != 0], diffT3[diffT3 != 0])/dim(data)[1]/10
    data$T1[sameT1T2] = data$T1[sameT1T2] - shift
    data$T1[sameT1T3] = data$T1[sameT1T3] - shift
    data$T2[sameT2T3] = data$T2[sameT2T3] - shift
    mindata = min(data$T1, data$T2, data$T3)
    if(mindata < 0){
      data$T1 = data$T1 + abs(mindata)
      data$T2 = data$T2 + abs(mindata)
      data$T3 = data$T3 + abs(mindata)
    }
    index_same_value = paste(sort(unique(c(sameT1T2, sameT1T3, sameT2T3))), collapse = ', ')
    if(push_warnings){warning("When T1 and T2 are observed, T1 should be less than T2. (Similar for T1/T3 and T2/T3).\nThe following number of rows: ", index_same_value,
                              " will be shifted by a small amount. Use $data to fetch the actual data. ")}
  }

  org_data = data
  data = data[sort(data$T3, index.return = TRUE)$ix, ]
  return(list(data = data, min_diff = min_diff, org_data = org_data))
}
#' @export
DNH_change_Z_to_12 = function(points, uZ){
  a = 1/(uZ[2] - uZ[1])
  b = 1 - uZ[1]/(uZ[2] - uZ[1])
  return(round(points * a + b))
}
#' @export
DNH_boot_shift_data = function(boot_data, min_diff){
  shift = runif(dim(boot_data)[1], min_diff/10, min_diff)
  boot_data$T1 = boot_data$T1 + shift
  boot_data$T2 = boot_data$T2 + shift
  boot_data$T3 = boot_data$T3 + shift
  return(boot_data)
}

# useful functions
#' @export
expit = function(x){1/(1 + exp(-x))}
#' @export
DNH_leave = function(){set.seed(Sys.time())}
#' @export
DNH_make_plot_main = function(i, case = NULL){
  if(is.null(case)){
    if(i == 1){main = expression(Delta[Z %->% Y])}
    if(i == 2){main = expression(Delta[Z %->% N[2] %->% Y])}
    if(i == 3){main = expression(Delta[Z %->% N[1] %->% Y])}
    if(i == 4){main = expression(Delta[Z %->% N[1] %->% N[2] %->% Y])}
  }else{
    if(i == 1){main = bquote('Case'~.(case)*':'~Delta[Z %->% Y])}
    if(i == 2){main = bquote('Case'~.(case)*':'~Delta[Z %->% N[2] %->% Y])}
    if(i == 3){main = bquote('Case'~.(case)*':'~Delta[Z %->% N[1] %->% Y])}
    if(i == 4){main = bquote('Case'~.(case)*':'~Delta[Z %->% N[1] %->% N[2] %->% Y])}
  }
  return(main)
}
#' @export
DNH_getcolor = function(gamma){
  if(gamma == 0){
    return("#000000")
  }else{
    gamma255 = floor(abs(gamma) * 255)
    gammacolor = toupper(as.hexmode(gamma255))
    if(gamma <= 0){
      return(paste("#0000", gammacolor, sep = ''))
    }else{
      return(paste("#", gammacolor, "0000", sep = ''))
    }
  }
}
#' @export
DNH_list2names = function(zvec){
  list_name = rep(NA, length(zvec))
  for(i in 1:length(zvec)){
    list_name[i] = paste("L", paste(zvec[[i]][[1]], collapse = ''), paste(zvec[[i]][[2]], collapse = ''), sep = '_')
  }
  return(list_name)
}
#' @export
DNH_create_case = function(gamma, effect_name, time){
  case = as.data.frame(matrix(0, nrow = length(time), ncol = length(gamma) + 1))
  colnames(case) = c('time', paste("gamma = ", gamma, sep = ''))
  case[, 1] = time
  case = list(case, case, case, case)
  names(case) = effect_name
  return(case)
}
#' @export
DNH_rep_row = function(x, n){
  return(matrix(rep(x, each = n), nrow = n))
}
#' @export
DNH_revcumsum = function(vec){
  return(rev(cumsum(rev(vec))))
}
#' @export
DNH_sort_mat = function(mat){
  mat_na = is.na(mat)
  if(sum(mat_na) == 0){
    mat = apply(mat, 2, sort)
  }else{
    mat = mat[rowSums(mat_na) == 0, ]
    mat = apply(mat, 2, sort)
  }
  return(mat)
}
#' @export
DNH_nan2zero = function(value){
  value[is.nan(value)] = 0
  return(value)
}
#' @export
DNH_inf2zero = function(value){
  value[is.infinite(value)] = 0
  return(value)
}
#' @export
DNH_get_position = function(x, y){
  # DNH_get_position returns vector z.
  # z_i := min_j(x_i<=y_j)
  # z_i := length(y) + 1, if x_i>y_j for all j
  # y is a ordered sequence. x may or may not be ordered. x and y can have different lengths.
  if(is.unsorted(y)){
    stop('The second argument must be sorted.')
  }
  z = approx(y, 1:length(y), x, yleft = 1, yright = length(y) + 1, f = 1, method = 'constant')$y
  return(z)
}
#' @export
DNH_get_sick_number = function(data, time_pt = NULL){
  if(is.null(dim(data))){
    data = t(as.matrix(data))
  }

  n_sample = dim(data)[1]
  n_col = dim(data)[2]
  if(n_col < 2){
    stop("number of data columns must >= 2")
  }

  if(n_col == 2){
    covariates = matrix(1, n_sample, 1)
  }else{
    covariates = as.matrix(data[, 3:n_col])
  }

  data = as.matrix(data[, 1:2], ncol = 2)
  if(!is.null(time_pt)){
    new_time = setdiff(time_pt, data[, 2])
    if(length(new_time) > 0){
      aug_data = matrix(0, nrow = dim(data)[1] + length(new_time), ncol = 2)
      aug_data[1:dim(data)[1], ] = data
      aug_data[dim(data)[1] + (1 : length(new_time)), ] = cbind(new_time, new_time)

      aug_cov = matrix(0, nrow = dim(data)[1] + length(new_time), ncol = dim(covariates)[2])
      aug_cov[1:dim(data)[1], ] = covariates
      aug_cov[dim(data)[1] + (1 : length(new_time)), ] = 0

      sort_id = sort(aug_data[, 2], index.return = TRUE)$ix
      data = aug_data[sort_id, ]
      covariates = as.matrix(aug_cov[sort_id, ], ncol = dim(aug_cov)[2])
    }
  }

  obs_time = matrix(data[, 2])
  time_data = data[, 1:2]
  if(dim(time_data)[1] == 1){time_data = t(time_data)}
  rank_time_2 = rank(time_data[, 2], ties.method = 'min')
  tmp_right = as.matrix(apply(covariates, 2, function(x) sum(x) - c(0, cumsum(x))[-(length(x) + 1)]))
  tmp_right = tmp_right[rank_time_2, ]

  sort_time_1 = sort(time_data[, 1], index.return = TRUE)
  sorted_time_df_1 = time_data[sort_time_1$ix, ]
  important_index = DNH_get_position(time_data[, 2], sort_time_1$x)

  sorted_covariates = as.matrix(covariates[sort_time_1$ix, ])
  tmp_left = apply(sorted_covariates, 2, function(x) sum(x) - c(0, cumsum(x)))
  tmp_left = tmp_left[important_index, ]

  if(!is.null(time_pt)){
    index_return = approx(data[, 2], 1:dim(data)[1], unique(time_pt), ties = 'max')$y
    if(is.null(dim(tmp_right))){
      return((tmp_right - tmp_left)[index_return])
    }else{
      return((tmp_right - tmp_left)[index_return, ])
    }
  }else{
    return(tmp_right - tmp_left)
  }
}

# functions help us get estimators
#' @export
DNH_get_Ybar = function(data, weights){
  dataZ1_index = (data$Z == 1)

  ## This part gives you YbarZ1_check (Z2) which is the same as YbarZ1 (Z2).
  # m = dim(data)[1]
  # dataZ1 = data[dataZ1_index, ]
  # dataZ2 = data[!dataZ1_index, ]
  # all_time = data$T3
  # YbarZ1_check = 0 * all_time
  # YbarZ2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   for(j in 1:m){
  #     if(j <= dim(dataZ1)[1]){YbarZ1_check[i] = YbarZ1_check[i] + (dataZ1$T3[j] >= all_time[i])}
  #     if(j <= dim(dataZ2)[1]){YbarZ2_check[i] = YbarZ2_check[i] + (dataZ2$T3[j] >= all_time[i])}
  #   }
  # }
  YbarZ1 = data.frame(time = data$T3, value = DNH_revcumsum(dataZ1_index * weights))
  YbarZ2 = data.frame(time = data$T3, value = DNH_revcumsum((!dataZ1_index) * weights))
  if(sum(diff(data$T3) == 0) > 0){
    YbarZ1$value = approx(YbarZ1$time, YbarZ1$value, data$T3, ties = 'max')$y
    YbarZ2$value = approx(YbarZ2$time, YbarZ2$value, data$T3, ties = 'max')$y
  }
  # YbarZ1$value - YbarZ1_check
  # YbarZ2$value - YbarZ2_check

  Ybar = list(YbarZ1 = YbarZ1, YbarZ2 = YbarZ2)
  return(Ybar)
}
#' @export
DNH_get_Ybar_n1.. = function(data, weights){
  dataZ1 = data[data$Z == 1, ]
  dataZ2 = data[data$Z == 2, ]

  weightsZ1 = weights[data$Z == 1]
  weightsZ2 = weights[data$Z == 2]

  ## n1 = 0
  ## This part gives you Ybar0.Z1_check (Z2) which is the same as Ybar0.Z1 (Z2).
  # all_time = data$T3
  # Ybar0.Z1_check = 0 * all_time
  # Ybar0.Z2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   for(j in 1:dim(data)[1]){
  #     if(j <= dim(dataZ1)[1]){Ybar0.Z1_check[i] = Ybar0.Z1_check[i] + ((dataZ1$T3[j] >= all_time[i]) & !(dataZ1$T1[j] < all_time[i] & dataZ1$d1[j]))}
  #     if(j <= dim(dataZ2)[1]){Ybar0.Z2_check[i] = Ybar0.Z2_check[i] + ((dataZ2$T3[j] >= all_time[i]) & !(dataZ2$T1[j] < all_time[i] & dataZ2$d1[j]))}
  #   }
  # }
  tmp_time = sort(c(dataZ1$T1[dataZ1$d1 == 1], dataZ1$T3[dataZ1$d1 == 0]))
  Ybar0.Z1 = data.frame(time = data$T3, value = approx(x = tmp_time, y = DNH_revcumsum(weightsZ1), xout = data$T3, method = 'constant', rule = 2, yright = 0, f = 1, ties = 'max')$y)
  tmp_time = sort(c(dataZ2$T1[dataZ2$d1 == 1], dataZ2$T3[dataZ2$d1 == 0]))
  Ybar0.Z2 = data.frame(time = data$T3, value = approx(x = tmp_time, y = DNH_revcumsum(weightsZ2), xout = data$T3, method = 'constant', rule = 2, yright = 0, f = 1, ties = 'max')$y)
  # Ybar0.Z1$value - Ybar0.Z1_check
  # Ybar0.Z2$value - Ybar0.Z2_check

  ## n1 = 1
  ## This part gives you Ybar1.Z1_check (Z2) which is the same as Ybar1.Z1 (Z2).
  # all_time = data$T3
  # Ybar1.Z1_check = 0 * all_time
  # Ybar1.Z2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   for(j in 1:dim(data)[1]){
  #     if(j <= dim(dataZ1)[1]){Ybar1.Z1_check[i] = Ybar1.Z1_check[i] + ((dataZ1$T3[j] >= all_time[i]) & (dataZ1$T1[j] < all_time[i] & dataZ1$d1[j]))}
  #     if(j <= dim(dataZ2)[1]){Ybar1.Z2_check[i] = Ybar1.Z2_check[i] + ((dataZ2$T3[j] >= all_time[i]) & (dataZ2$T1[j] < all_time[i] & dataZ2$d1[j]))}
  #   }
  # }
  if(sum(dataZ1$d1 == 1) == 0){
    Ybar1.Z1 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar1.Z1_tmp = DNH_get_sick_number(cbind(dataZ1[dataZ1$d1 == 1, c(1, 5)], weightsZ1[dataZ1$d1 == 1]), data$T3)
    if(length(Ybar1.Z1_tmp) == length(data$T3)){
      Ybar1.Z1 = data.frame(time = data$T3, value = Ybar1.Z1_tmp)
    }else{
      Ybar1.Z1 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar1.Z1_tmp, data$T3, method = 'constant')$y)
    }
  }
  if(sum(dataZ2$d1 == 1) == 0){
    Ybar1.Z2 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar1.Z2_tmp = DNH_get_sick_number(cbind(dataZ2[dataZ2$d1 == 1, c(1, 5)], weightsZ2[dataZ2$d1 == 1]), data$T3)
    if(length(Ybar1.Z2_tmp) == length(data$T3)){
      Ybar1.Z2 = data.frame(time = data$T3, value = Ybar1.Z2_tmp)
    }else{
      Ybar1.Z2 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar1.Z2_tmp, data$T3, method = 'constant')$y)
    }
  }

  # Ybar1.Z1$value - Ybar1.Z1_check
  # Ybar1.Z2$value - Ybar1.Z2_check

  Ybarn1.. = list(Ybar0.Z1 = Ybar0.Z1, Ybar0.Z2 = Ybar0.Z2, Ybar1.Z1 = Ybar1.Z1, Ybar1.Z2 = Ybar1.Z2)
  return(Ybarn1..)
}
#' @export
DNH_get_Ybar_..n2 = function(data, weights){
  dataZ1 = data[data$Z == 1, ]
  dataZ2 = data[data$Z == 2, ]

  weightsZ1 = weights[data$Z == 1]
  weightsZ2 = weights[data$Z == 2]

  ## n2 = 0
  ## This part gives you Ybar.0Z1_check (Z2) which is the same as Ybar.0Z1 (Z2).
  # all_time = data$T3
  # Ybar.0Z1_check = 0 * all_time
  # Ybar.0Z2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   for(j in 1:dim(data)[1]){
  #     if(j <= dim(dataZ1)[1]){Ybar.0Z1_check[i] = Ybar.0Z1_check[i] + ((dataZ1$T3[j] >= all_time[i]) & !(dataZ1$T2[j] < all_time[i] & dataZ1$d2[j]))}
  #     if(j <= dim(dataZ2)[1]){Ybar.0Z2_check[i] = Ybar.0Z2_check[i] + ((dataZ2$T3[j] >= all_time[i]) & !(dataZ2$T2[j] < all_time[i] & dataZ2$d2[j]))}
  #   }
  # }
  Ybar.0Z1 = data.frame(time = data$T3, value = approx(x = sort(dataZ1$T2), y = DNH_revcumsum(weightsZ1), xout = data$T3, method = 'constant', rule = 2, yright = 0, f = 1, ties = 'max')$y)
  Ybar.0Z2 = data.frame(time = data$T3, value = approx(x = sort(dataZ2$T2), y = DNH_revcumsum(weightsZ2), xout = data$T3, method = 'constant', rule = 2, yright = 0, f = 1, ties = 'max')$y)
  # Ybar.0Z1$value - Ybar.0Z1_check
  # Ybar.0Z2$value - Ybar.0Z2_check

  ## n2 = 1
  ## This part gives you Ybar.1Z1_check (Z2) which is the same as Ybar.1Z1 (Z2).
  # all_time = data$T3
  # Ybar.1Z1_check = 0 * all_time
  # Ybar.1Z2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   for(j in 1:dim(data)[1]){
  #     if(j <= dim(dataZ1)[1]){Ybar.1Z1_check[i] = Ybar.1Z1_check[i] + ((dataZ1$T3[j] >= all_time[i]) & (dataZ1$T2[j] < all_time[i] & dataZ1$d2[j]))}
  #     if(j <= dim(dataZ2)[1]){Ybar.1Z2_check[i] = Ybar.1Z2_check[i] + ((dataZ2$T3[j] >= all_time[i]) & (dataZ2$T2[j] < all_time[i] & dataZ2$d2[j]))}
  #   }
  # }
  if(sum(dataZ1$d2 == 1) == 0){
    Ybar.1Z1 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar.1Z1_tmp = DNH_get_sick_number(cbind(dataZ1[dataZ1$d2 == 1, c(3, 5)], weightsZ1[dataZ1$d2 == 1]), data$T3)
    if(length(Ybar.1Z1_tmp) == length(data$T3)){
      Ybar.1Z1 = data.frame(time = data$T3, value = Ybar.1Z1_tmp)
    }else{
      Ybar.1Z1 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar.1Z1_tmp, data$T3, method = 'constant')$y)
    }
  }

  if(sum(dataZ2$d2 == 1) == 0){
    Ybar.1Z2 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar.1Z2_tmp = DNH_get_sick_number(cbind(dataZ2[dataZ2$d2 == 1, c(3, 5)], weightsZ2[dataZ2$d2 == 1]), data$T3)
    if(length(Ybar.1Z2_tmp) == length(data$T3)){
      Ybar.1Z2 = data.frame(time = data$T3, value = Ybar.1Z2_tmp)
    }else{
      Ybar.1Z2 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar.1Z2_tmp, data$T3, method = 'constant')$y)
    }
  }

  # Ybar.1Z1$value - Ybar.1Z1_check
  # Ybar.1Z2$value - Ybar.1Z2_check

  Ybar..n2 = list(Ybar.0Z1 = Ybar.0Z1, Ybar.0Z2 = Ybar.0Z2, Ybar.1Z1 = Ybar.1Z1, Ybar.1Z2 = Ybar.1Z2)
  return(Ybar..n2)
}
#' @export
DNH_get_Ybar_n1n2 = function(data, Ybar, Ybarn1.., Ybar..n2, weights){
  dataZ1 = data[data$Z == 1, ]
  dataZ2 = data[data$Z == 2, ]

  weightsZ1 = weights[data$Z == 1]
  weightsZ2 = weights[data$Z == 2]

  ## n1 = 1, n2 = 1
  # all_time = data$T3
  # Ybar11Z1_check = 0 * all_time
  # Ybar11Z2_check = 0 * all_time
  # for(i in 1:length(all_time)){
  #   print(i)
  #   for(j in 1:dim(data)[1]){
  #     if(j <= dim(dataZ1)[1]){Ybar11Z1_check[i] = Ybar11Z1_check[i] + ((dataZ1$T3[j] >= all_time[i]) & (dataZ1$T1[j] < all_time[i] & dataZ1$d1[j]) & (dataZ1$T2[j] < all_time[i] & dataZ1$d2[j]))}
  #     if(j <= dim(dataZ2)[1]){Ybar11Z2_check[i] = Ybar11Z2_check[i] + ((dataZ2$T3[j] >= all_time[i]) & (dataZ2$T1[j] < all_time[i] & dataZ2$d1[j]) & (dataZ2$T2[j] < all_time[i] & dataZ2$d2[j]))}
  #   }
  # }
  if(sum(dataZ1$d1 == 1 & dataZ1$d2 == 1) == 0){
    Ybar11Z1 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar11Z1_tmp = DNH_get_sick_number(cbind(dataZ1[dataZ1$d1 == 1 & dataZ1$d2 == 1, c(3, 5)], weightsZ1[dataZ1$d1 == 1 & dataZ1$d2 == 1]), data$T3)
    if(length(Ybar11Z1_tmp) == length(data$T3)){
      Ybar11Z1 = data.frame(time = data$T3, value = Ybar11Z1_tmp)
    }else{
      Ybar11Z1 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar11Z1_tmp, data$T3, method = 'constant')$y)
    }
  }

  if(sum(dataZ2$d1 == 1 & dataZ2$d2 == 1) == 0){
    Ybar11Z2 = data.frame(time = data$T3, value = 0)
  }else{
    Ybar11Z2_tmp = DNH_get_sick_number(cbind(dataZ2[dataZ2$d1 == 1 & dataZ2$d2 == 1, c(3, 5)], weightsZ2[dataZ2$d1 == 1 & dataZ2$d2 == 1]), data$T3)
    if(length(Ybar11Z2_tmp) == length(data$T3)){
      Ybar11Z2 = data.frame(time = data$T3, value = Ybar11Z2_tmp)
    }else{
      Ybar11Z2 = data.frame(time = data$T3, value = approx(unique(data$T3), Ybar11Z2_tmp, data$T3, method = 'constant')$y)
    }
  }

  # Ybar11Z1$value - Ybar11Z1_check
  # Ybar11Z2$value - Ybar11Z2_check
  ## n1 = 0, n2 = 1
  Ybar01Z1 = data.frame(time = data$T3, value = Ybar..n2$Ybar.1Z1$value - Ybar11Z1$value)
  Ybar01Z2 = data.frame(time = data$T3, value = Ybar..n2$Ybar.1Z2$value - Ybar11Z2$value)

  ## n1 = 1, n2 = 0
  Ybar10Z1 = data.frame(time = data$T3, value = Ybarn1..$Ybar1.Z1$value - Ybar11Z1$value)
  Ybar10Z2 = data.frame(time = data$T3, value = Ybarn1..$Ybar1.Z2$value - Ybar11Z2$value)

  ## n1 = 0, n2 = 0
  Ybar00Z1 = data.frame(time = data$T3, value = Ybar$YbarZ1$value - Ybarn1..$Ybar1.Z1$value - Ybar..n2$Ybar.1Z1$value + Ybar11Z1$value)
  Ybar00Z2 = data.frame(time = data$T3, value = Ybar$YbarZ2$value - Ybarn1..$Ybar1.Z2$value - Ybar..n2$Ybar.1Z2$value + Ybar11Z2$value)

  Ybarn1n2 = list(Ybar00Z1 = Ybar00Z1, Ybar00Z2 = Ybar00Z2, Ybar01Z1 = Ybar01Z1, Ybar01Z2 = Ybar01Z2, Ybar10Z1 = Ybar10Z1, Ybar10Z2 = Ybar10Z2, Ybar11Z1 = Ybar11Z1, Ybar11Z2 = Ybar11Z2)
}
#' @export
DNH_get_w_n1n2 = function(Ybar, Ybarn1n2){
  time = Ybar$YbarZ1$time

  ## n1 = 0, n2 = 0
  w00 = list(Z1 = NULL, Z2 = NULL)
  w00$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar00Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1n2$Ybar00Z1$time]))
  w00$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar00Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1n2$Ybar00Z2$time]))

  ## n1 = 0, n2 = 1
  w01 = list(Z1 = NULL, Z2 = NULL)
  w01$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar01Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1n2$Ybar01Z1$time]))
  w01$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar01Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1n2$Ybar01Z2$time]))

  ## n1 = 1, n2 = 0
  w10 = list(Z1 = NULL, Z2 = NULL)
  w10$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar10Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1n2$Ybar10Z1$time]))
  w10$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar10Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1n2$Ybar10Z2$time]))

  ## n1 = 1, n2 = 1
  w11 = list(Z1 = NULL, Z2 = NULL)
  w11$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar11Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1n2$Ybar11Z1$time]))
  w11$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1n2$Ybar11Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1n2$Ybar11Z2$time]))

  wn1n2 = list(w00 = w00, w01 = w01, w10 = w10, w11 = w11)
  return(wn1n2)
  # plot(Ybar$YbarZ1$time, approx(c(0, w00$Z1[[1]]), c(1, w00$Z1[[2]]), Ybar$YbarZ1$time, rule = 2)$y + approx(w01$Z1[[1]], w01$Z1[[2]], Ybar$YbarZ1$time, rule = 2)$y + approx(w10$Z1[[1]], w10$Z1[[2]], Ybar$YbarZ1$time, rule = 2)$y + approx(w11$Z1[[1]], w11$Z1[[2]], Ybar$YbarZ1$time, rule = 2)$y, ylab = '')
  # points(Ybar$YbarZ2$time, approx(c(0, w00$Z2[[1]]), c(1, w00$Z2[[2]]), Ybar$YbarZ2$time, rule = 2)$y + approx(w01$Z2[[1]], w01$Z2[[2]], Ybar$YbarZ2$time, rule = 2)$y + approx(w10$Z2[[1]], w10$Z2[[2]], Ybar$YbarZ2$time, rule = 2)$y + approx(w11$Z2[[1]], w11$Z2[[2]], Ybar$YbarZ2$time, rule = 2)$y, ylab = '', col = 'red')
}
#' @export
DNH_get_w_n1.. = function(wn1n2, Ybar, Ybarn1..){
  time = Ybar$YbarZ1$time

  ## n1 = 0
  w0. = list(Z1 = NULL, Z2 = NULL)
  w0.$Z1 = data.frame(time = time, value = wn1n2$w00$Z1$value + wn1n2$w01$Z1$value)
  w0.$Z2 = data.frame(time = time, value = wn1n2$w00$Z2$value + wn1n2$w01$Z2$value)
  # w0.$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1..$Ybar0.Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1..$Ybar0.Z1$time]))
  # w0.$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1..$Ybar0.Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1..$Ybar0.Z2$time]))

  ## n1 = 1
  w1. = list(Z1 = NULL, Z2 = NULL)
  w1.$Z1 = data.frame(time = time, value = wn1n2$w10$Z1$value + wn1n2$w11$Z1$value)
  w1.$Z2 = data.frame(time = time, value = wn1n2$w10$Z2$value + wn1n2$w11$Z2$value)
  # w1.$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybarn1..$Ybar1.Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybarn1..$Ybar1.Z1$time]))
  # w1.$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybarn1..$Ybar1.Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybarn1..$Ybar1.Z2$time]))

  wn1.. = list(w0. = w0., w1. = w1.)
  return(wn1..)
  # plot(Ybar$YbarZ1$time, approx(c(0, w0.$Z1[[1]]), c(1, w0.$Z1[[2]]), Ybar$YbarZ1$time, rule = 2)$y + approx(c(0, w1.$Z1[[1]]), c(0, w1.$Z1[[2]]), Ybar$YbarZ1$time, rule = 2)$y, ylab = '')
  # points(Ybar$YbarZ2$time, approx(c(0, w0.$Z2[[1]]), c(1, w0.$Z2[[2]]), Ybar$YbarZ2$time, rule = 2)$y + approx(c(0, w1.$Z2[[1]]), c(0, w1.$Z2[[2]]), Ybar$YbarZ2$time, rule = 2)$y, ylab = '', col = 'red')
}
#' @export
DNH_get_w_..n2 = function(wn1n2, Ybar, Ybar..n2){
  time = Ybar$YbarZ1$time

  ## n2 = 0
  w.0 = list(Z1 = NULL, Z2 = NULL)
  w.0$Z1 = data.frame(time = time, value = wn1n2$w00$Z1$value + wn1n2$w10$Z1$value)
  w.0$Z2 = data.frame(time = time, value = wn1n2$w00$Z2$value + wn1n2$w10$Z2$value)
  # w.0$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybar..n2$Ybar.0Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybar..n2$Ybar.0Z1$time]))
  # w.0$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybar..n2$Ybar.0Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybar..n2$Ybar.0Z2$time]))

  ## n2 = 1
  w.1 = list(Z1 = NULL, Z2 = NULL)
  w.1$Z1 = data.frame(time = time, value = wn1n2$w01$Z1$value + wn1n2$w11$Z1$value)
  w.1$Z2 = data.frame(time = time, value = wn1n2$w01$Z2$value + wn1n2$w11$Z2$value)
  # w.1$Z1 = data.frame(time = time, value = DNH_nan2zero(Ybar..n2$Ybar.1Z1$value / Ybar$YbarZ1$value[Ybar$YbarZ1$time %in% Ybar..n2$Ybar.1Z1$time]))
  # w.1$Z2 = data.frame(time = time, value = DNH_nan2zero(Ybar..n2$Ybar.1Z2$value / Ybar$YbarZ2$value[Ybar$YbarZ2$time %in% Ybar..n2$Ybar.1Z2$time]))

  w..n2 = list(w.0 = w.0, w.1 = w.1)
  return(w..n2)

  # plot(Ybar$YbarZ1$time, approx(c(0, w.0$Z1[[1]]), c(1, w.0$Z1[[2]]), Ybar$YbarZ1$time, rule = 2)$y + approx(c(0, w.1$Z1[[1]]), c(0, w.1$Z1[[2]]), Ybar$YbarZ1$time, rule = 2)$y, ylab = '')
  # points(Ybar$YbarZ2$time, approx(c(0, w.0$Z2[[1]]), c(1, w.0$Z2[[2]]), Ybar$YbarZ2$time, rule = 2)$y + approx(c(0, w.1$Z2[[1]]), c(0, w.1$Z2[[2]]), Ybar$YbarZ2$time, rule = 2)$y, ylab = '', col = 'red')
}
#' @export
DNH_get_wcond_n1n2 = function(Ybarn1n2, Ybarn1..){
  ## n1 = 0, n2 = 0
  wcond00 = list(Z1 = NULL, Z2 = NULL)
  wcond00$Z1 = data.frame(time = Ybarn1n2$Ybar00Z1$time, value = DNH_nan2zero(Ybarn1n2$Ybar00Z1$value / Ybarn1..$Ybar0.Z1$value[Ybarn1..$Ybar0.Z1$time %in% Ybarn1n2$Ybar00Z1$time]))
  wcond00$Z2 = data.frame(time = Ybarn1n2$Ybar00Z2$time, value = DNH_nan2zero(Ybarn1n2$Ybar00Z2$value / Ybarn1..$Ybar0.Z2$value[Ybarn1..$Ybar0.Z2$time %in% Ybarn1n2$Ybar00Z2$time]))

  ## n1 = 0, n2 = 1
  wcond01 = list(Z1 = NULL, Z2 = NULL)
  wcond01$Z1 = data.frame(time = Ybarn1n2$Ybar10Z1$time, value = DNH_nan2zero(Ybarn1n2$Ybar10Z1$value / Ybarn1..$Ybar1.Z1$value[Ybarn1..$Ybar1.Z1$time %in% Ybarn1n2$Ybar10Z1$time]))
  wcond01$Z2 = data.frame(time = Ybarn1n2$Ybar10Z2$time, value = DNH_nan2zero(Ybarn1n2$Ybar10Z2$value / Ybarn1..$Ybar1.Z2$value[Ybarn1..$Ybar1.Z2$time %in% Ybarn1n2$Ybar10Z2$time]))

  ## n1 = 1, n2 = 0
  wcond10 = list(Z1 = NULL, Z2 = NULL)
  wcond10$Z1 = data.frame(time = Ybarn1n2$Ybar01Z1$time, value = DNH_nan2zero(Ybarn1n2$Ybar01Z1$value / Ybarn1..$Ybar0.Z1$value[Ybarn1..$Ybar0.Z1$time %in% Ybarn1n2$Ybar01Z1$time]))
  wcond10$Z2 = data.frame(time = Ybarn1n2$Ybar01Z2$time, value = DNH_nan2zero(Ybarn1n2$Ybar01Z2$value / Ybarn1..$Ybar0.Z2$value[Ybarn1..$Ybar0.Z2$time %in% Ybarn1n2$Ybar01Z2$time]))

  ## n1 = 1, n2 = 1
  wcond11 = list(Z1 = NULL, Z2 = NULL)
  wcond11$Z1 = data.frame(time = Ybarn1n2$Ybar11Z1$time, value = DNH_nan2zero(Ybarn1n2$Ybar11Z1$value / Ybarn1..$Ybar1.Z1$value[Ybarn1..$Ybar1.Z1$time %in% Ybarn1n2$Ybar11Z1$time]))
  wcond11$Z2 = data.frame(time = Ybarn1n2$Ybar11Z2$time, value = DNH_nan2zero(Ybarn1n2$Ybar11Z2$value / Ybarn1..$Ybar1.Z2$value[Ybarn1..$Ybar1.Z2$time %in% Ybarn1n2$Ybar11Z2$time]))

  wcondn1n2 = list(wcond00 = wcond00, wcond01 = wcond01, wcond10 = wcond10, wcond11 = wcond11)
  return(wcondn1n2)

  # plot(Ybarn1..$Ybar0.Z1$time, approx(c(0, wcond00$Z1[[1]]), c(1, wcond00$Z1[[2]]), Ybarn1..$Ybar0.Z1$time, rule = 2)$y + approx(c(0, wcond10$Z1[[1]]), c(0, wcond10$Z1[[2]]), Ybarn1..$Ybar0.Z1$time, rule = 2)$y, ylab = '', ylim = c(0.9, 1.1))
  # plot(Ybarn1..$Ybar0.Z2$time, approx(c(0, wcond00$Z2[[1]]), c(1, wcond00$Z2[[2]]), Ybarn1..$Ybar0.Z2$time, rule = 2)$y + approx(c(0, wcond10$Z2[[1]]), c(0, wcond10$Z2[[2]]), Ybarn1..$Ybar0.Z2$time, rule = 2)$y, ylab = '', ylim = c(0.9, 1.1))
  # plot(Ybarn1..$Ybar1.Z1$time, approx(c(0, wcond01$Z1[[1]]), c(1, wcond01$Z1[[2]]), Ybarn1..$Ybar1.Z1$time, rule = 2)$y + approx(c(0, wcond11$Z1[[1]]), c(0, wcond11$Z1[[2]]), Ybarn1..$Ybar1.Z1$time, rule = 2)$y, ylab = '', ylim = c(0.9, 1.1))
  # plot(Ybarn1..$Ybar1.Z1$time, approx(c(0, wcond01$Z2[[1]]), c(1, wcond01$Z2[[2]]), Ybarn1..$Ybar1.Z1$time, rule = 2)$y + approx(c(0, wcond11$Z2[[1]]), c(0, wcond11$Z2[[2]]), Ybarn1..$Ybar1.Z1$time, rule = 2)$y, ylab = '', ylim = c(0.9, 1.1))
}

# estimators and counterfactual hazard
#' @export
DNH_get_weights = function(dLn1n2, iPSEweights = TRUE){
  time = dLn1n2$dL00$Z1$time
  len = length(time)
  if(iPSEweights){
    ## n1 = n2 = 0
    weights00 = list(Z1 = NULL, Z2 = NULL)
    weights00$Z1 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL00$Z1$value[-len]))))
    weights00$Z2 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL00$Z2$value[-len]))))

    ## n1 = 0, n2 = 1
    weights01 = list(Z1 = NULL, Z2 = NULL)
    weights01$Z1 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL01$Z1$value[-len]))))
    weights01$Z2 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL01$Z2$value[-len]))))

    ## n1 = 1, n2 = 0
    weights10 = list(Z1 = NULL, Z2 = NULL)
    weights10$Z1 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL10$Z1$value[-len]))))
    weights10$Z2 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL10$Z2$value[-len]))))

    ## n1 = 1, n2 = 1
    weights11 = list(Z1 = NULL, Z2 = NULL)
    weights11$Z1 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL11$Z1$value[-len]))))
    weights11$Z2 = data.frame(time = time, value = c(1, exp(-cumsum(dLn1n2$dL11$Z2$value[-len]))))
  }else{
    tmp = data.frame(time = time, value = rep(1, length(time)))
    weights00 = list(Z1 = tmp, Z2 = tmp)
    weights01 = list(Z1 = tmp, Z2 = tmp)
    weights10 = list(Z1 = tmp, Z2 = tmp)
    weights11 = list(Z1 = tmp, Z2 = tmp)
  }

  weights = list(weights00 = weights00, weights01 = weights01, weights10 = weights10, weights11 = weights11)
  return(weights)
}
#' @export
DNH_get_dN3bar_n1n2 = function(data, weights){
  m = dim(data)[1]
  data = data[sort(data$T3, index.return = TRUE)$ix, ]
  dataZ1 = data[data$Z == 1, ]
  dataZ2 = data[data$Z == 2, ]

  ## n1 = 0, n2 = 0
  boolean00 = data$d1 == 0 & data$d2 == 0 & data$d3 == 1
  dN3bar00Z1 = rep(0, m)
  dN3bar00Z2 = rep(0, m)
  dN3bar00Z1[boolean00 & data$Z == 1] = weights[boolean00 & data$Z == 1]
  dN3bar00Z2[boolean00 & data$Z == 2] = weights[boolean00 & data$Z == 2]
  dN3bar00Z1 = data.frame(time = data$T3, value = dN3bar00Z1)
  dN3bar00Z2 = data.frame(time = data$T3, value = dN3bar00Z2)

  ## n1 = 0, n2 = 1
  boolean01 = data$d1 == 0 & data$d2 == 1 & data$d3 == 1
  dN3bar01Z1 = rep(0, m)
  dN3bar01Z2 = rep(0, m)
  dN3bar01Z1[boolean01 & data$Z == 1] = weights[boolean01 & data$Z == 1]
  dN3bar01Z2[boolean01 & data$Z == 2] = weights[boolean01 & data$Z == 2]
  dN3bar01Z1 = data.frame(time = data$T3, value = dN3bar01Z1)
  dN3bar01Z2 = data.frame(time = data$T3, value = dN3bar01Z2)

  ## n1 = 1, n2 = 0
  boolean10 = data$d1 == 1 & data$d2 == 0 & data$d3 == 1
  dN3bar10Z1 = rep(0, m)
  dN3bar10Z2 = rep(0, m)
  dN3bar10Z1[boolean10 & data$Z == 1] = weights[boolean10 & data$Z == 1]
  dN3bar10Z2[boolean10 & data$Z == 2] = weights[boolean10 & data$Z == 2]
  dN3bar10Z1 = data.frame(time = data$T3, value = dN3bar10Z1)
  dN3bar10Z2 = data.frame(time = data$T3, value = dN3bar10Z2)

  ## n1 = 1, n2 = 1
  boolean11 = data$d1 == 1 & data$d2 == 1 & data$d3 == 1
  dN3bar11Z1 = rep(0, m)
  dN3bar11Z2 = rep(0, m)
  dN3bar11Z1[boolean11 & data$Z == 1] = weights[boolean11 & data$Z == 1]
  dN3bar11Z2[boolean11 & data$Z == 2] = weights[boolean11 & data$Z == 2]
  dN3bar11Z1 = data.frame(time = data$T3, value = dN3bar11Z1)
  dN3bar11Z2 = data.frame(time = data$T3, value = dN3bar11Z2)

  dN3barn1n2 = list(dN3bar00Z1 = dN3bar00Z1, dN3bar00Z2 = dN3bar00Z2, dN3bar01Z1 = dN3bar01Z1, dN3bar01Z2 = dN3bar01Z2,
                    dN3bar10Z1 = dN3bar10Z1, dN3bar10Z2 = dN3bar10Z2, dN3bar11Z1 = dN3bar11Z1, dN3bar11Z2 = dN3bar11Z2)
  return(dN3barn1n2)
}
#' @export
DNH_get_wa = function(iPSEweights, components, wbn1n2, zvec){
  za = zvec[1]
  time = wbn1n2$wb00$time
  if(iPSEweights){weights = components$iPSEweights}
  if(!iPSEweights){weights = components$NoiPSEweights}
  denominator = wbn1n2$wb00$value * weights$weights00[[za]]$value + wbn1n2$wb01$value * weights$weights01[[za]]$value + wbn1n2$wb10$value * weights$weights10[[za]]$value + wbn1n2$wb11$value * weights$weights11[[za]]$value
  wa00 = data.frame(time = time, value = weights$weights00[[za]]$value/denominator)
  wa01 = data.frame(time = time, value = weights$weights01[[za]]$value/denominator)
  wa10 = data.frame(time = time, value = weights$weights10[[za]]$value/denominator)
  wa11 = data.frame(time = time, value = weights$weights11[[za]]$value/denominator)

  zero_deno = (denominator == 0)
  if(sum(zero_deno)){
    sum_weights = (weights$weights00[[za]]$value + weights$weights01[[za]]$value + weights$weights10[[za]]$value + weights$weights11[[za]]$value)[zero_deno]
    wa00$value[zero_deno] = weights$weights00[[za]]$value[zero_deno]/sum_weights
    wa01$value[zero_deno] = weights$weights01[[za]]$value[zero_deno]/sum_weights
    wa10$value[zero_deno] = weights$weights10[[za]]$value[zero_deno]/sum_weights
    wa11$value[zero_deno] = weights$weights11[[za]]$value[zero_deno]/sum_weights
  }
  wan1n2 = list(wa00 = wa00, wa01 = wa01, wa10 = wa10, wa11 = wa11)
  return(wan1n2)
}
#' @export
DNH_get_wb = function(wcondn1n2, wn1.., zvec){
  zb = zvec[2]; zc = zvec[3]; zd = zvec[4]

  if(FALSE){
    wb00 = data.frame(time = wn1..$w0.[[zc]]$time, value = wcondn1n2$wcond00[[zb]]$value * wn1..$w0.[[zc]]$value)
    wb01 = data.frame(time = wn1..$w0.[[zc]]$time, value = wcondn1n2$wcond10[[zb]]$value * wn1..$w0.[[zc]]$value)
    wb10 = data.frame(time = wn1..$w1.[[zc]]$time, value = wcondn1n2$wcond01[[zb]]$value * wn1..$w1.[[zc]]$value)
    wb11 = data.frame(time = wn1..$w1.[[zc]]$time, value = wcondn1n2$wcond11[[zb]]$value * wn1..$w1.[[zc]]$value)
  }else{
    wb00 = data.frame(time = wn1..$w0.[[zc]]$time,
                      value = (wcondn1n2$wcond00[[zb]]$value * wn1..$w0.[[zd]]$value + wcondn1n2$wcond01[[zb]]$value * wn1..$w1.[[zd]]$value) * wn1..$w0.[[zc]]$value)
    wb01 = data.frame(time = wn1..$w0.[[zc]]$time,
                      value = (wcondn1n2$wcond10[[zb]]$value * wn1..$w0.[[zd]]$value + wcondn1n2$wcond11[[zb]]$value * wn1..$w1.[[zd]]$value) * wn1..$w0.[[zc]]$value)
    wb10 = data.frame(time = wn1..$w1.[[zc]]$time,
                      value = (wcondn1n2$wcond00[[zb]]$value * wn1..$w0.[[zd]]$value + wcondn1n2$wcond01[[zb]]$value * wn1..$w1.[[zd]]$value) * wn1..$w1.[[zc]]$value)
    wb11 = data.frame(time = wn1..$w1.[[zc]]$time,
                      value = (wcondn1n2$wcond10[[zb]]$value * wn1..$w1.[[zd]]$value + wcondn1n2$wcond11[[zb]]$value * wn1..$w1.[[zd]]$value) * wn1..$w1.[[zc]]$value)
  }

  wbn1n2 = list(wb00 = wb00, wb01 = wb01, wb10 = wb10, wb11 = wb11)
  return(wbn1n2)
}
#' @export
DNH_get_dL_n1n2 = function(dN3barn1n2, Ybarn1n2){
  dL00 = list(Z1 = NULL, Z2 = NULL)
  dL00$Z1 = data.frame(time = dN3barn1n2$dN3bar00Z1$time, value = DNH_nan2zero(dN3barn1n2$dN3bar00Z1$value/Ybarn1n2$Ybar00Z1$value))
  dL00$Z2 = data.frame(time = dN3barn1n2$dN3bar00Z2$time, value = DNH_nan2zero(dN3barn1n2$dN3bar00Z2$value/Ybarn1n2$Ybar00Z2$value))
  dL01 = list(Z1 = NULL, Z2 = NULL)
  dL01$Z1 = data.frame(time = dN3barn1n2$dN3bar01Z1$time, value = DNH_nan2zero(dN3barn1n2$dN3bar01Z1$value/Ybarn1n2$Ybar01Z1$value))
  dL01$Z2 = data.frame(time = dN3barn1n2$dN3bar01Z2$time, value = DNH_nan2zero(dN3barn1n2$dN3bar01Z2$value/Ybarn1n2$Ybar01Z2$value))
  dL10 = list(Z1 = NULL, Z2 = NULL)
  dL10$Z1 = data.frame(time = dN3barn1n2$dN3bar10Z1$time, value = DNH_nan2zero(dN3barn1n2$dN3bar10Z1$value/Ybarn1n2$Ybar10Z1$value))
  dL10$Z2 = data.frame(time = dN3barn1n2$dN3bar10Z2$time, value = DNH_nan2zero(dN3barn1n2$dN3bar10Z2$value/Ybarn1n2$Ybar10Z2$value))
  dL11 = list(Z1 = NULL, Z2 = NULL)
  dL11$Z1 = data.frame(time = dN3barn1n2$dN3bar11Z1$time, value = DNH_nan2zero(dN3barn1n2$dN3bar11Z1$value/Ybarn1n2$Ybar11Z1$value))
  dL11$Z2 = data.frame(time = dN3barn1n2$dN3bar11Z2$time, value = DNH_nan2zero(dN3barn1n2$dN3bar11Z2$value/Ybarn1n2$Ybar11Z2$value))

  dLn1n2 = list(dL00 = dL00, dL01 = dL01, dL10 = dL10, dL11 = dL11)
  return(dLn1n2)
}
#' @export
DNH_get_components = function(data, weights){
  m = dim(data)[1]

  ## Y and mu related
  Ybar = DNH_get_Ybar(data, weights)
  Ybarn1.. = DNH_get_Ybar_n1..(data, weights)
  Ybar..n2 = DNH_get_Ybar_..n2(data, weights)
  Ybarn1n2 = DNH_get_Ybar_n1n2(data, Ybar, Ybarn1.., Ybar..n2, weights)

  ## w related
  wn1n2 = DNH_get_w_n1n2(Ybar, Ybarn1n2)
  wn1.. = DNH_get_w_n1..(wn1n2, Ybar, Ybarn1..)
  w..n2 = DNH_get_w_..n2(wn1n2, Ybar, Ybar..n2)
  wcondn1n2 = DNH_get_wcond_n1n2(Ybarn1n2, Ybarn1..)

  ## dL related
  dN3barn1n2 = DNH_get_dN3bar_n1n2(data, weights)
  dLn1n2 = DNH_get_dL_n1n2(dN3barn1n2, Ybarn1n2)

  ## weights related
  iPSEweights = DNH_get_weights(dLn1n2, iPSEweights = TRUE)
  NoiPSEweights = DNH_get_weights(dLn1n2, iPSEweights = FALSE)

  ## index
  observed = Ybar$YbarZ1$time %in% data$T3[data$d3 == 1]
  components = list(observed = observed, wn1.. = wn1.., w..n2 = w..n2, wcondn1n2 = wcondn1n2, iPSEweights = iPSEweights, dLn1n2 = dLn1n2, NoiPSEweights = NoiPSEweights,
                    Ybar = Ybar, dN3barn1n2 = dN3barn1n2, Ybarn1n2 = Ybarn1n2, Ybarn1.. = Ybarn1.., Ybar..n2 = Ybar..n2, wn1n2 = wn1n2)
  return(components)
}
#' @export
DNH_get_sen_components = function(logit_dLn1n2, components, bU, gZ, gN1, gN2){
  sen_components = components
  sen_components$dLn1n2$dL00$Z1$value = expit(logit_dLn1n2$dL00$Z1$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 0 + gN1 * 0 + gN2 * 0))
  sen_components$dLn1n2$dL00$Z2$value = expit(logit_dLn1n2$dL00$Z2$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 1 + gN1 * 0 + gN2 * 0))
  sen_components$dLn1n2$dL01$Z1$value = expit(logit_dLn1n2$dL01$Z1$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 0 + gN1 * 0 + gN2 * 1))
  sen_components$dLn1n2$dL01$Z2$value = expit(logit_dLn1n2$dL01$Z2$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 1 + gN1 * 0 + gN2 * 1))
  sen_components$dLn1n2$dL10$Z1$value = expit(logit_dLn1n2$dL10$Z1$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 0 + gN1 * 1 + gN2 * 0))
  sen_components$dLn1n2$dL10$Z2$value = expit(logit_dLn1n2$dL10$Z2$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 1 + gN1 * 1 + gN2 * 0))
  sen_components$dLn1n2$dL11$Z1$value = expit(logit_dLn1n2$dL11$Z1$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 0 + gN1 * 1 + gN2 * 1))
  sen_components$dLn1n2$dL11$Z2$value = expit(logit_dLn1n2$dL11$Z2$value * sqrt(1 + 0.35 * bU^2) - bU * (gZ * 1 + gN1 * 1 + gN2 * 1))
  sen_components$iPSEweights = DNH_get_weights(sen_components$dLn1n2)

  return(sen_components)
}
#' @export
DNH_get_counterfactual_hazard = function(wan1n2, wbn1n2, dLn1n2, zvec){
  za = zvec[1]

  Lambda00 = cumsum(wan1n2$wa00$value * wbn1n2$wb00$value * dLn1n2$dL00[[za]]$value)
  Lambda01 = cumsum(wan1n2$wa01$value * wbn1n2$wb01$value * dLn1n2$dL01[[za]]$value)
  Lambda10 = cumsum(wan1n2$wa10$value * wbn1n2$wb10$value * dLn1n2$dL10[[za]]$value)
  Lambda11 = cumsum(wan1n2$wa11$value * wbn1n2$wb11$value * dLn1n2$dL11[[za]]$value)

  Lambda = data.frame(time = wbn1n2$wb00$time, value = Lambda00 + Lambda01 + Lambda10 + Lambda11)
  return(Lambda)
}

# effect
#' @export
DNH_get_effect = function(components, zvec, iPSEweights){
  if(is.null(zvec)){
    zveclist = list(c(2, 1, 1, 1), c(1, 1, 1, 1),
                    c(2, 2, 1, 1), c(2, 1, 1, 1),
                    c(2, 2, 2, 1), c(2, 2, 1, 1),
                    c(2, 2, 2, 2), c(2, 2, 2, 1))
    effect_tmp = vector(mode = 'list', length = 8)
    for(i in 1:8){
      zvec = zveclist[[i]]
      wbn1n2 = DNH_get_wb(components$wcondn1n2, components$wn1.., zvec)
      wan1n2 = DNH_get_wa(iPSEweights, components, wbn1n2, zvec)
      effect_tmp[[i]] = DNH_get_counterfactual_hazard(wan1n2, wbn1n2, components$dLn1n2, zvec)
    }

    ## get effect
    effect = list(iPSEZY = NULL, iPSEZ2Y = NULL, iPSEZ1Y = NULL, iPSEZ12Y = NULL)
    effect[[1]] = data.frame(time = effect_tmp[[1]]$time, effect = effect_tmp[[1]]$value - effect_tmp[[2]]$value)
    effect[[2]] = data.frame(time = effect_tmp[[1]]$time, effect = effect_tmp[[3]]$value - effect_tmp[[4]]$value)
    effect[[3]] = data.frame(time = effect_tmp[[1]]$time, effect = effect_tmp[[5]]$value - effect_tmp[[6]]$value)
    effect[[4]] = data.frame(time = effect_tmp[[1]]$time, effect = effect_tmp[[7]]$value - effect_tmp[[8]]$value)
  }else{
    m = length(components$Ybar$YbarZ1$time)
    effect = vector(mode = 'list', length = length(zvec))
    names(effect) = DNH_list2names(zvec)
    for(i in 1:length(zvec)){
      effect[[i]] = data.frame(time = components$Ybar$YbarZ1$time, effect = rep(0, m))
      for(j in 1:2){
        zvec_now = DNH_change_Z_to_12(zvec[[i]][[j]], c(0, 1))
        wbn1n2 = DNH_get_wb(components$wcondn1n2, components$wn1.., zvec_now)
        if(iPSEweights){fn1n2 = DNH_get_f_n1n2(components$iPSEweights, wbn1n2, zvec_now)}
        if(!iPSEweights){fn1n2 = DNH_get_f_n1n2(components$NoiPSEweights, wbn1n2, zvec_now)}
        bigWn1n2 = DNH_get_bigW_n1n2_ver2(fn1n2)
        effect[[i]]$effect = effect[[i]]$effect + (-1)^(j - 1) * DNH_get_counterfactual_hazard(bigWn1n2, components$dLn1n2, zvec_now)$value
      }
    }
  }

  return(effect)
}
#' @export
DNH_get_sen_ana = function(data, components, zvec, iPSEweights, effect_name){
  # effect_name = names(effect)
  effect = list(case1 = NULL, case2 = NULL, case3 = NULL, case4 = NULL, case5 = NULL)
  logit_dLn1n2 = components$dLn1n2
  logit_dLn1n2$dL00$Z1$value = pracma::logit(components$dLn1n2$dL00$Z1$value)
  logit_dLn1n2$dL00$Z2$value = pracma::logit(components$dLn1n2$dL00$Z2$value)
  logit_dLn1n2$dL01$Z1$value = pracma::logit(components$dLn1n2$dL01$Z1$value)
  logit_dLn1n2$dL01$Z2$value = pracma::logit(components$dLn1n2$dL01$Z2$value)
  logit_dLn1n2$dL10$Z1$value = pracma::logit(components$dLn1n2$dL10$Z1$value)
  logit_dLn1n2$dL10$Z2$value = pracma::logit(components$dLn1n2$dL10$Z2$value)
  logit_dLn1n2$dL11$Z1$value = pracma::logit(components$dLn1n2$dL11$Z1$value)
  logit_dLn1n2$dL11$Z2$value = pracma::logit(components$dLn1n2$dL11$Z2$value)

  ## case 1
  gamma = seq(-1, 1, length.out = 11)
  effect$case1 = DNH_create_case(gamma, time = components$Ybar$YbarZ1$time, effect_name)
  betaU = gamma
  gammaZ = abs(gamma)
  gammaN1 = abs(gamma)
  gammaN2 = abs(gamma)
  for(i in 1:length(gamma)){
    bU = betaU[i]; gZ = gammaZ[i]; gN1 = gammaN1[i]; gN2 = gammaN2[i]
    sen_components = DNH_get_sen_components(logit_dLn1n2, components, bU, gZ, gN1, gN2)
    tmp = DNH_get_effect(sen_components, zvec, iPSEweights)
    for(j in 1:length(effect_name)){effect$case1[[j]][[i + 1]] = tmp[[j]]$effect}
  }

  ## case 2
  gamma = seq(0, 1, length.out = 11)
  effect$case2 = DNH_create_case(gamma, time = components$Ybar$YbarZ1$time, effect_name)
  betaU = log(1.5)
  gammaZ = abs(gamma)
  gammaN1 = abs(gamma)
  gammaN2 = abs(gamma)
  for(i in 1:length(gamma)){
    bU = betaU; gZ = gammaZ[i]; gN1 = gammaN1[i]; gN2 = gammaN2[i]
    sen_components = DNH_get_sen_components(logit_dLn1n2, components, bU, gZ, gN1, gN2)
    tmp = DNH_get_effect(sen_components, zvec, iPSEweights)
    for(j in 1:length(effect_name)){effect$case2[[j]][[i + 1]] = tmp[[j]]$effect}
  }

  ## case 3
  gamma = seq(0, 1, length.out = 11)
  effect$case3 = DNH_create_case(gamma, time = components$Ybar$YbarZ1$time, effect_name)
  betaU = -log(1.5)
  gammaZ = abs(gamma)
  gammaN1 = abs(gamma)
  gammaN2 = abs(gamma)
  for(i in 1:length(gamma)){
    bU = betaU; gZ = gammaZ[i]; gN1 = gammaN1[i]; gN2 = gammaN2[i]
    sen_components = DNH_get_sen_components(logit_dLn1n2, components, bU, gZ, gN1, gN2)
    tmp = DNH_get_effect(sen_components, zvec, iPSEweights)
    for(j in 1:length(effect_name)){effect$case3[[j]][[i + 1]] = tmp[[j]]$effect}
  }

  ## case 4
  gamma = seq(-1, 1, length.out = 11)
  effect$case4 = DNH_create_case(gamma, time = components$Ybar$YbarZ1$time, effect_name)
  betaU = gamma
  gammaZ = 0.2
  gammaN1 = 0.2
  gammaN2 = 0.2
  for(i in 1:length(gamma)){
    bU = betaU[i]; gZ = gammaZ; gN1 = gammaN1; gN2 = gammaN2
    sen_components = DNH_get_sen_components(logit_dLn1n2, components, bU, gZ, gN1, gN2)
    tmp = DNH_get_effect(sen_components, zvec, iPSEweights)
    for(j in 1:length(effect_name)){effect$case4[[j]][[i + 1]] = tmp[[j]]$effect}
  }

  ## case 5
  gamma = seq(-1, 1, length.out = 11)
  effect$case5 = DNH_create_case(gamma, time = components$Ybar$YbarZ1$time, effect_name)
  betaU = 0
  gammaZ = gamma
  gammaN1 = gamma
  gammaN2 = gamma
  for(i in 1:length(gamma)){
    bU = betaU; gZ = gammaZ[i]; gN1 = gammaN1[i]; gN2 = gammaN2[i]
    sen_components = DNH_get_sen_components(logit_dLn1n2, components, bU, gZ, gN1, gN2)
    tmp = DNH_get_effect(sen_components, zvec, iPSEweights)
    for(j in 1:length(effect_name)){effect$case5[[j]][[i + 1]] = tmp[[j]]$effect}
  }
  return(effect)
}
#' @export
DNH_PM = function(effect, PM){
  time_len = length(effect[[1]]$time)
  effect_len = length(effect)
  total_effect = effect[[1]]$time * 0
  if(PM == TRUE){for(i in 1:effect_len){total_effect = total_effect + effect[[i]]$effect}}
  if(PM == 'abs'){for(i in 1:effect_len){total_effect = total_effect + abs(effect[[i]]$effect)}}
  total_effect_integral = sum((total_effect[2:time_len] + total_effect[1:(time_len - 1)])/2 * diff(effect[[i]]$time))
  proportion_mediated = data.frame(rep(0, effect_len), row.names = names(effect))
  if(PM == TRUE){for(i in 1:effect_len){proportion_mediated[i, 1] = sum((effect[[i]]$effect[2:time_len] + effect[[i]]$effect[1:(time_len - 1)])/2 * diff(effect[[i]]$time))/total_effect_integral}}
  if(PM == 'abs'){for(i in 1:effect_len){proportion_mediated[i, 1] = sum((abs(effect[[i]]$effect[2:time_len]) + abs(effect[[i]]$effect[1:(time_len - 1)]))/2 * diff(effect[[i]]$time))/total_effect_integral}}
  return(proportion_mediated)
}

# main function
#' @export
DNH_plot_result = function(effect, zvec, plot_unit, match_ylim){
  xlab = ifelse(floor(plot_unit) == 365, 'Time (years)', 'time')
  plot_cut = floor(length(effect[[1]]$time) * 0.99)
  time_min = effect[[1]]$time[1]/plot_unit
  time = effect[[1]]$time[plot_cut]/plot_unit

  if(match_ylim){
    ylimnow = c(Inf, -Inf)
    for(i in 1:length(effect)){
      ylimnow[1] = min(ylimnow[1], min(effect[[i]]$boot_lower[1:plot_cut]))
      ylimnow[2] = max(ylimnow[2], max(effect[[i]]$boot_upper[1:plot_cut]))
    }
  }
  lwd = c(1.5, 2, 2, 2)
  col_effect = rep('black', 4)
  col_boot = rep('blue', 4)
  # col_effect = c('deeppink3', 'royalblue1', 'purple3', 'orange3')
  # col_boot = c('deeppink', 'steelblue1', 'purple', 'sandybrown')
  for(i in 1:length(effect)){
    if(is.null(zvec)){
      main = DNH_make_plot_main(i)
    }else{
      main = paste('iPSE, (', paste(unlist(zvec[[i]][1]), collapse = ', '), ') - (', paste(unlist(zvec[[i]][2]), collapse = ', '), ')', sep = '')
    }
    if(!match_ylim){ylimnow = c(min(effect[[i]]$boot_lower[1:plot_cut]), max(effect[[i]]$boot_upper[1:plot_cut]))}
    plot(effect[[i]]$time/plot_unit, effect[[i]]$effect, xlim = c(time_min, time), ylim = ylimnow, col = col_effect[i], type = 's', xlab = xlab, ylab = 'effect', main = main, cex.main = 1.5, cex.lab = 1.2, lwd = lwd[i])
    lines(effect[[i]]$time/plot_unit, effect[[i]]$boot_lower, col = col_boot[i], lwd = lwd[i])
    lines(effect[[i]]$time/plot_unit, effect[[i]]$boot_upper, col = col_boot[i], lwd = lwd[i])
    abline(h = 0, col = 'grey')
  }
}
#' @export
DNH_save_plot_result = function(effect, zvec, plot_unit, plot_name, folder_name, match_ylim){
  # plot_unit = 365.25
  # folder_name = '~/Huang_DNH/RealData/plot_result/test'
  # plot_name = 'tmp'
  if(is.null(folder_name)){folder_name = '~'}
  xlab = ifelse(floor(plot_unit) == 365, 'Time (years)', 'time')

  plot_cut = floor(length(effect[[1]]$time) * 0.99)
  time_min = effect[[1]]$time[1]/plot_unit
  time = effect[[1]]$time[plot_cut]/plot_unit

  if(match_ylim){
    ylimnow = c(Inf, -Inf)
    for(i in 1:length(effect)){
      ylimnow[1] = min(ylimnow[1], min(effect[[i]]$boot_lower[1:plot_cut]))
      ylimnow[2] = max(ylimnow[2], max(effect[[i]]$boot_upper[1:plot_cut]))
    }
  }
  lwd = c(1.5, 2, 2, 2)
  col_effect = rep('black', 4)
  col_boot = rep('blue', 4)
  # col_effect = c('deeppink3', 'royalblue1', 'purple3', 'orange3')
  # col_boot = c('deeppink', 'steelblue1', 'purple', 'sandybrown')
  for(i in 1:length(effect)){
    if(i == 1){png(paste(folder_name, '/', plot_name, '_ZY.png', sep = ''), pointsize = 16)}
    if(i == 2){png(paste(folder_name, '/', plot_name, '_Z2Y.png', sep = ''), pointsize = 16)}
    if(i == 3){png(paste(folder_name, '/', plot_name, '_Z1Y.png', sep = ''), pointsize = 16)}
    if(i == 4){png(paste(folder_name, '/', plot_name, '_Z12Y.png', sep = ''), pointsize = 16)}
    if(is.null(zvec)){
      main = DNH_make_plot_main(i)
    }else{
      main = paste('iPSE, (', paste(unlist(zvec[[i]][1]), collapse = ', '), ') - (', paste(unlist(zvec[[i]][2]), collapse = ', '), ')', sep = '')
    }
    if(!match_ylim){ylimnow = c(min(effect[[i]]$boot_lower[1:plot_cut]), max(effect[[i]]$boot_upper[1:plot_cut]))}
    plot(effect[[i]]$time/plot_unit, effect[[i]]$effect, xlim = c(time_min, time), ylim = ylimnow, col = col_effect[i], type = 's', xlab = xlab, ylab = 'effect', main = main, cex.main = 1.5, cex.lab = 1.2, lwd = lwd[i])
    lines(effect[[i]]$time/plot_unit, effect[[i]]$boot_lower, col = col_boot[i], lwd = lwd[i])
    lines(effect[[i]]$time/plot_unit, effect[[i]]$boot_upper, col = col_boot[i], lwd = lwd[i])
    abline(h = 0, col = 'grey')
    dev.off()
  }
}
#' @export
DNH_plot_sen_ana = function(sen_ana, zvec, plot_unit, sensitivity_analysis_match_ylim){
  xlab = ifelse(floor(plot_unit) == 365, 'Time (years)', 'time')
  time = sen_ana$case1$iPSEZY$time/plot_unit
  plot_cut = floor(length(time) * 0.99)
  time_min = time[1]

  ylim_min = matrix(0, nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  ylim_max = matrix(0, nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  for(i in 1:length(sen_ana)){
    for(j in 1:length(sen_ana[[1]])){
      ylim_min[i, j] = min(sen_ana[[i]][[j]][, 2:(dim(sen_ana[[i]][[j]])[2])])
      ylim_max[i, j] = max(sen_ana[[i]][[j]][, 2:(dim(sen_ana[[i]][[j]])[2])])
    }
  }
  if(sensitivity_analysis_match_ylim == TRUE){
    ylim_min = matrix(min(ylim_min), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
    ylim_max = matrix(max(ylim_max), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  }else if(sensitivity_analysis_match_ylim == 'case'){
    ylim_min = matrix(apply(ylim_min, 1, min), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
    ylim_max = matrix(apply(ylim_max, 1, max), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  }else if(sensitivity_analysis_match_ylim == 'effect'){
    ylim_min = DNH_rep_row(apply(ylim_min, 2, min), length(sen_ana))
    ylim_max = DNH_rep_row(apply(ylim_max, 2, max), length(sen_ana))
  }
  ylim_min = ylim_min * 1.5; ylim_max = ylim_max * 1.5;

  for(i in 1:length(sen_ana)){
    for(j in 1:length(sen_ana[[1]])){
      if(is.null(zvec)){
        main = DNH_make_plot_main(j, case = i)
      }else{
        main = paste('Case ', i, ': iPSE, (', paste(unlist(zvec[[i]][1]), collapse = ', '), ') - (', paste(unlist(zvec[[i]][2]), collapse = ', '), ')', sep = '')
      }
      now_gamma = rep(0, 11); now_col = rep(NA, 11)
      now_gamma[1] = -1 + (i == 2) + (i == 3)
      now_col[1] = DNH_getcolor(now_gamma[1])
      lwd = 1 + (now_gamma[1] == 0)
      plot(time, sen_ana[[i]][[j]][, 2], ylim = c(ylim_min[i, j], ylim_max[i, j] * 1.5), col = now_col[1], type = 's', xlab = xlab, ylab = 'effect', main = main, cex.main = 1.5, cex.lab = 1.2, lwd = lwd)
      for(k in 3:dim(sen_ana[[i]][[j]])[2]){
        now_gamma[k - 1] = 0.9/(dim(sen_ana[[i]][[j]])[2] - 3) * (k - 3) + 0.1
        if(i == 1 || i == 4 || i == 5){now_gamma[k - 1] = 2 * (now_gamma[k - 1] - 0.1) - 0.8}
        now_col[k - 1] = DNH_getcolor(now_gamma[k - 1])
        lwd = 1 + (now_gamma[k - 1] == 0)
        lines(time, sen_ana[[i]][[j]][, k], col = now_col[k - 1], lwd = lwd)
      }
      if(mean(as.matrix(sen_ana[[i]][[j]][, 2:dim(sen_ana[[i]][[j]])[2]]), na.rm = TRUE) < (ylim_min[i, j] + ylim_max[i, j]) / 2){
        legend("topleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2, bty = 'n')
      }else{
        legend("bottomleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2, bty = 'n')
      }

      abline(h = 0, col = 'grey')
    }
  }
}
#' @export
DNH_save_plot_sen_ana = function(sen_ana, zvec, plot_unit, plot_name, folder_name, sensitivity_analysis_match_ylim){
  # plot_unit = 365.25
  # folder_name = '~/Huang_DNH/RealData/plot_result/test'
  # plot_name = 'tmp'
  xlab = ifelse(floor(plot_unit) == 365, 'Time (years)', 'time')
  time = sen_ana$case1$iPSEZY$time/plot_unit
  plot_cut = floor(length(time) * 0.99)
  time_min = time[1]

  ylim_min = matrix(0, nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  ylim_max = matrix(0, nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  for(i in 1:length(sen_ana)){
    for(j in 1:length(sen_ana[[1]])){
      ylim_min[i, j] = min(sen_ana[[i]][[j]][, 2:(dim(sen_ana[[i]][[j]])[2])])
      ylim_max[i, j] = max(sen_ana[[i]][[j]][, 2:(dim(sen_ana[[i]][[j]])[2])])
    }
  }
  if(sensitivity_analysis_match_ylim == TRUE){
    ylim_min = matrix(min(ylim_min), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
    ylim_max = matrix(max(ylim_max), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  }else if(sensitivity_analysis_match_ylim == 'case'){
    ylim_min = matrix(apply(ylim_min, 1, min), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
    ylim_max = matrix(apply(ylim_max, 1, max), nrow = length(sen_ana), ncol = length(sen_ana[[1]]))
  }else if(sensitivity_analysis_match_ylim == 'effect'){
    ylim_min = DNH_rep_row(apply(ylim_min, 2, min), length(sen_ana))
    ylim_max = DNH_rep_row(apply(ylim_max, 2, max), length(sen_ana))
  }
  ylim_min = ylim_min * 1.5; ylim_max = ylim_max * 1.5;

  for(i in 1:length(sen_ana)){
    for(j in 1:length(sen_ana[[1]])){
      if(j == 1){png(paste(folder_name, '/', 'Case', i, '_', plot_name, '_ZY.png', sep = ''), pointsize = 16)}
      if(j == 2){png(paste(folder_name, '/', 'Case', i, '_', plot_name, '_Z2Y.png', sep = ''), pointsize = 16)}
      if(j == 3){png(paste(folder_name, '/', 'Case', i, '_', plot_name, '_Z1Y.png', sep = ''), pointsize = 16)}
      if(j == 4){png(paste(folder_name, '/', 'Case', i, '_', plot_name, '_Z12Y.png', sep = ''), pointsize = 16)}
      if(is.null(zvec)){
        main = DNH_make_plot_main(j, case = i)
      }else{
        main = paste('Case ', i, ': iPSE, (', paste(unlist(zvec[[i]][1]), collapse = ', '), ') - (', paste(unlist(zvec[[i]][2]), collapse = ', '), ')', sep = '')
      }
      now_gamma = rep(0, 11); now_col = rep(NA, 11)
      now_gamma[1] = -1 + (i == 2) + (i == 3)
      now_col[1] = DNH_getcolor(now_gamma[1])
      lwd = 1 + (now_gamma[1] == 0)
      plot(time, sen_ana[[i]][[j]][, 2], ylim = c(ylim_min[i, j], ylim_max[i, j]), col = now_col[1], type = 's', xlab = xlab, ylab = 'effect', main = main, cex.main = 1.5, cex.lab = 1.2, lwd = lwd)
      for(k in 3:dim(sen_ana[[i]][[j]])[2]){
        now_gamma[k - 1] = 0.9/(dim(sen_ana[[i]][[j]])[2] - 3) * (k - 3) + 0.1
        if(i == 1 || i == 4 || i == 5){now_gamma[k - 1] = 2 * (now_gamma[k - 1] - 0.1) - 0.8}
        now_col[k - 1] = DNH_getcolor(now_gamma[k - 1])
        lwd = 1 + (now_gamma[k - 1] == 0)
        lines(time, sen_ana[[i]][[j]][, k], col = now_col[k - 1], lwd = lwd)
      }
      if(mean(as.matrix(sen_ana[[i]][[j]][, 2:dim(sen_ana[[i]][[j]])[2]]), na.rm = TRUE) < (ylim_min[i, j] + ylim_max[i, j]) / 2){
        legend("topleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2, bty = 'n')
      }else{
        legend("bottomleft", legend = now_gamma, col = now_col, lty = 1, lwd = 2, title = expression(gamma), ncol = 2, bty = 'n')
      }
      abline(h = 0, col = 'grey')
      dev.off()
    }
  }
}

#' @export
DNH_Med_TTEM = function(components){
  time = components$Ybar$YbarZ1$time[components$observed]

  ## T1
  resultM1 = list(DE = NULL, IE = NULL)
  # (1, 1)
  prob_n0 = components$wn1..$w0.$Z2$value
  dLambda_n0 = DNH_nan2zero((components$dN3barn1n2$dN3bar00Z2$value + components$dN3barn1n2$dN3bar01Z2$value)/components$Ybarn1..$Ybar0.Z2$value)
  prob_n1 = components$wn1..$w1.$Z2$value
  dLambda_n1 = DNH_nan2zero((components$dN3barn1n2$dN3bar10Z2$value + components$dN3barn1n2$dN3bar11Z2$value)/components$Ybarn1..$Ybar1.Z2$value)
  tmp1 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)

  # (1, 0)
  prob_n0 = components$wn1..$w0.$Z1$value
  prob_n1 = components$wn1..$w1.$Z1$value
  tmp2 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)

  # (0, 0)
  dLambda_n0 = DNH_nan2zero((components$dN3barn1n2$dN3bar00Z1$value + components$dN3barn1n2$dN3bar01Z1$value)/components$Ybarn1..$Ybar0.Z1$value)
  dLambda_n1 = DNH_nan2zero((components$dN3barn1n2$dN3bar10Z1$value + components$dN3barn1n2$dN3bar11Z1$value)/components$Ybarn1..$Ybar1.Z1$value)
  tmp3 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)
  resultM1$DE = data.frame(time = time, value = (tmp2 - tmp3)[components$observed])
  resultM1$IE = data.frame(time = time, value = (tmp1 - tmp2)[components$observed])

  ## T2
  resultM2 = list(DE = NULL, IE = NULL)
  # (1, 1)
  prob_n0 = components$w..n2$w.0$Z2$value
  dLambda_n0 = DNH_nan2zero((components$dN3barn1n2$dN3bar00Z2$value + components$dN3barn1n2$dN3bar10Z2$value)/components$Ybar..n2$Ybar.0Z2$value)
  prob_n1 = components$w..n2$w.1$Z2$value
  dLambda_n1 = DNH_nan2zero((components$dN3barn1n2$dN3bar01Z2$value + components$dN3barn1n2$dN3bar11Z2$value)/components$Ybar..n2$Ybar.1Z2$value)
  tmp1 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)

  # (1, 0)
  prob_n0 = components$w..n2$w.0$Z1$value
  prob_n1 = components$w..n2$w.1$Z1$value
  tmp2 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)

  # (0, 0)
  dLambda_n0 = DNH_nan2zero((components$dN3barn1n2$dN3bar00Z1$value + components$dN3barn1n2$dN3bar10Z1$value)/components$Ybar..n2$Ybar.0Z1$value)
  dLambda_n1 = DNH_nan2zero((components$dN3barn1n2$dN3bar01Z1$value + components$dN3barn1n2$dN3bar11Z1$value)/components$Ybar..n2$Ybar.1Z1$value)
  tmp3 = cumsum(prob_n0 * dLambda_n0 + prob_n1 * dLambda_n1)
  resultM2$DE = data.frame(time = time, value = (tmp2 - tmp3)[components$observed])
  resultM2$IE = data.frame(time = time, value = (tmp1 - tmp2)[components$observed])

  result = list(M1 = resultM1, M2 = resultM2)
  return(result)
}
