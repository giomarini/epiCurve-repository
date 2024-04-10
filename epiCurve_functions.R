#### epi_fit #####
# This function fits the observed epidemiological curves providing an estimate for the three free parameters
# of the FOI function.
# It uses the below likelihood function to find the best combination maximizing the likelihood
epi_fit=function(observed_dataset,population_dataset,output_file_name){
  
  matrix_weekly_cases=observed_dataset
  population=population_dataset
  
  days=1:364
  
  columns_names=paste0("W",1:52)
  
  # creating the output dataframe
  output_model=data.frame(REGIONS=matrix_weekly_cases$REGION,
                          YEAR=matrix_weekly_cases$YEAR)
  output_model$POPULATION=NA
  output_model$OBSERVED_TOTAL=NA
  output_model$MODEL_TOTAL=NA
  output_model$MU=NA
  output_model$SD=NA
  output_model$C_NORM=NA
  
  pb <- txtProgressBar(min = 0, max = nrow(matrix_weekly_cases), style = 3)
  for(i in 1:nrow(matrix_weekly_cases)){
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
    
    weekly_observed_cases=as.numeric(matrix_weekly_cases[i,columns_names])
    selected_region=matrix_weekly_cases$REGION[i]
    N=as.numeric(population$POP[which(population$REGION==selected_region)])
    
    #looking for a starting point
    n_iter=1
    curr_lik=Inf
    set.seed(123)
    while(n_iter<5000){
      c_FOI_tmp=runif(n=1,min=0,max=0.001)
      mu_FOI_tmp=runif(n=1,min=160,max=300)
      sd_FOI_tmp=runif(n=1,min=1,max=50)
      lik_tmp=pois_likelihood_func(weekly_observed_cases,N,c(c_FOI_tmp,mu_FOI_tmp,sd_FOI_tmp))
      if(lik_tmp<curr_lik){
        parms_kept=c(c_FOI_tmp,mu_FOI_tmp,sd_FOI_tmp)
        curr_lik=lik_tmp
      } 
      n_iter=n_iter+1
    }
    
    # pois_likelihood_func(weekly_observed_cases,N,parms_kept)
    
    if(is.finite(curr_lik)){
      optimized=optim(par=parms_kept, fn=pois_likelihood_func, 
                      observed_dataset=weekly_observed_cases,pop=N)
      
      weekly_cases_model=c()
      for(week_start_index in seq(1,max(days),7)){
        c_FOI=optimized$par[1]
        mu_FOI=optimized$par[2]
        sd_FOI=optimized$par[3]
        FOI_model=c_FOI*dnorm(days, mean = mu_FOI, sd = sd_FOI)
        
        tmp=0
        for(giorno in week_start_index:(week_start_index+6))
          tmp=tmp+N*FOI_model[giorno]
        weekly_cases_model=c(weekly_cases_model,tmp)
      }
      weekly_cases_model=round(weekly_cases_model)
      
      output_model$POPULATION[i]=N
      
      output_model$OBSERVED_TOTAL[i]=sum(weekly_observed_cases,na.rm=T)
      output_model$MODEL_TOTAL[i]=sum(weekly_cases_model,na.rm=T)
      
      output_model$C_NORM[i]=optimized$par[1]
      output_model$MU[i]=round(optimized$par[2],2)
      output_model$SD[i]=round(optimized$par[3],2)
    }
  }
  write.csv(output_model,file=output_file_name,row.names=F)
}

#### pois_likelihood_func ####
# It is assumed cases follow a Poisson distribution
pois_likelihood_func=function(observed_dataset,pop,FOI_params){
  days=1:364
  
  c_FOI=FOI_params[1]
  mu_FOI=FOI_params[2]
  sd_FOI=FOI_params[3]
  FOI_to_test=c_FOI*dnorm(days, mean = mu_FOI, sd = sd_FOI)
  N=pop
  
  cum_log_prob=0
  dataset_index=1
  for(week_start_index in seq(1,max(days),7)){
    weekly_cases=0
    for(day in week_start_index:(week_start_index+6)){
      daily_cases=N*FOI_to_test[day]
      weekly_cases=weekly_cases+daily_cases
    }
    weekly_cases=round(weekly_cases)
    prob_week=dpois(x=observed_dataset[dataset_index],lambda=weekly_cases,log=T)
    cum_log_prob=cum_log_prob+prob_week
    dataset_index=dataset_index+1
  }
  return(-cum_log_prob)
}

##### create_dataset ####
# This function creates the aggregated dataset (with weekly cases) used in the epi_fit function 
# starting from the raw case_based dataset.
# minimum_cases is the minimum number of cases for an epidemiological curve to be modeled
create_dataset=function(minimum_cases,cases_df){
  # creating YEAR and WEEK variables
  cases_df$YEAR=year(cases_df$ONSET_DATE)
  cases_df$WEEK=week(cases_df$ONSET_DATE)
  
  total_region_year=cases_df %>%
    group_by(REGION,YEAR) %>%
    summarise(TOTAL_CASES=n())
  selected=which(total_region_year$TOTAL_CASES>=minimum_cases)
  
  matrix_weekly_cases=matrix(NA,nrow=length(selected),ncol=2+52)
  for(i in 1:length(selected)){
    selected_row=selected[i]
    selected_region=total_region_year$REGION[selected_row]
    selected_year=total_region_year$YEAR[selected_row]
    cum_weekly_tmp=rep(0,52)
    for(j in 1:52){
      cum_weekly_tmp[j]=length(which(cases_df$REGION==selected_region &
                                       cases_df$YEAR==selected_year &
                                       cases_df$WEEK==j))
    }
    matrix_weekly_cases[i,]=c(selected_region,selected_year,cum_weekly_tmp)
  }
  matrix_weekly_cases=as.data.frame(matrix_weekly_cases)
  for(j in 3:ncol(matrix_weekly_cases))
    matrix_weekly_cases[,j]=as.numeric(matrix_weekly_cases[,j])
  names(matrix_weekly_cases)=c("REGION","YEAR",paste0("W",c(1:52)))
  
  return(matrix_weekly_cases)
}
