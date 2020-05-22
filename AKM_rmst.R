
# --- RMST Using Adjusted KM ---
# Time is the time to event
# Status is 0 if censored, 1 if event
# Group should be a factor variable
# Weights can be obtained separately, ie through logistic models
# Tau is a user-specified truncation point. 
# If not specified, the default will be the minimum of the each groups' last event time 

akm_rmst <- function(time, status, group, weight=NULL, tau=NULL, alpha=.05, 
                     xaxismin=0, xaxismax=max(time)){
  if(sum(time<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=0)>0){print("Error: weights must be greater than 0.")
    }else{
      if(sum(status!=0 & status!=1)>0){print("Error: status must be a vector of 0s and/or 1s.")
      }else{
        
        if(is.null(weight)){weight <- rep(1, length(time))}	
        data <- data.frame(time, status, group, weight)
        data <- data[!is.na(data$group) & !is.na(data$time),]
        data <- data[order(group),] 

        #--- If tau not specified, use minimum tau from all groups ---
        j=length(unique(data$group))
        
        if(is.null(tau)){
          taui = rep(999, j)
          for (i in (1:j)){
            groupval <- (levels(data$group)[i])
            dat_group <- data[which(data$group==(groupval)),]
            taui[i] <- max(dat_group$time[dat_group$status==1])
            }
           tau <- min(taui)
         }

        #--- Calculate AKM RMST in each group ---
        rmst <- rep(999, length(1:j))
        groupval <- rep(999, length(1:j))
        rmst_var <- rep(999, length(1:j))
        rmst_se <- rep(999, length(1:j))
        plot(NULL, xlim=c(xaxismin, xaxismax), ylim=c(0,1), xlab='Time',ylab='Adjusted Survival Probability')
        title(main='Adjusted Kaplan-Meier')

        for (i in 1:j){
          groupval[i] <- (levels(data$group)[i])
          dat_group <- data[which(data$group==(groupval[i])),]
          
          #--- AKM ---
          # Based on 'adjusted.KM' function from {IPWsurvival} package
          # Author: F. Le Borgne and Y. Foucher
          tj <- c(0,sort(unique(dat_group$time[dat_group$status==1])))
          dj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time==x & dat_group$status==1])})
          yj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time>=x])})
          st <- cumprod(1-(dj/yj))
          m <- sapply(tj, function(x){sum((dat_group$weight[dat_group$time>=x])^2)})
          mj <- ((yj^2)/m)
          #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
          ft <- data.frame(tj, yj, dj, st, i, mj)
          
          #--- RMST ---
          # Based on 'rmst1 function' from {survRM2} package
          # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
          rtime <- ft$tj<=tau
          tj_r <- sort(c(ft$tj[rtime],tau))
          st_r <- ft$st[rtime]
          yj_r <- ft$yj[rtime]
          dj_r <- ft$dj[rtime]
          time_diff <- diff(c(0, tj_r))
          areas <- time_diff * c(1, st_r)
          rmst[i] <- sum(areas)
          
          #--- Variance ---
          mj_r <- ft$mj[rtime]
          var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
          #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
          var_r <- c(var_r,0)
          rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
          rmst_se[i] <- sqrt(rmst_var[i])
          
          #--- Plot AKM ---
          lines(ft$tj, ft$st,type="s", col=(i+2), lwd=2)
        }
      }
    }
  }
  
  #--- Add legend and tau to plot ---
  abline(v=tau, col=1, lty=3, lwd=2)
  legend('bottomleft', paste("Group", groupval), lty=rep(1, j), lwd=rep(2, j), col=3:(j+2), 
         cex=.75, bty ="n", inset = c(0, 0))
  
  #--- Compare RMST between groups and compile output---
  results <- data.frame(groupval,rmst,rmst_var,rmst_se,tau)
  pwc <- ((j^2)-j)/2   #number of pairwise comparisons
  
  label_diff <- rep(999,(pwc))
  rmst_diff <- rep(999,(pwc))
  rmst_diff_se <- rep(999,(pwc))
  rmst_diff_low <- rep(999,(pwc))
  rmst_diff_upp <- rep(999,(pwc))
  rmst_diff_pval <- rep(999,(pwc))
  
  label_ratio <- rep(999,(pwc))
  rmst_logratio <- rep(999,(pwc))
  rmst_logratio_se <- rep(999,(pwc))
  rmst_ratio <- rep(999,(pwc))
  rmst_ratio_low <- rep(999,(pwc))
  rmst_ratio_upp <- rep(999,(pwc))
  rmst_logratio_pval <- rep(999,(pwc))
  
  output_diff <- data.frame(label_diff,rmst_diff,rmst_diff_se,rmst_diff_low,rmst_diff_upp,rmst_diff_pval)
  output_ratio <- data.frame(label_ratio,rmst_logratio,rmst_logratio_se,rmst_ratio,rmst_ratio_low,rmst_ratio_upp,rmst_logratio_pval)
  l <- 1
  
  for (i in 1:(j-1)){
    for (j in (i+1):j){
      # Based on 'rmst2 function' from {survRM2} package
      # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
      
      #--- RMST Difference ---
      output_diff[l,]$label_diff <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_diff[l,]$rmst_diff <- (results[j,]$rmst - results[i,]$rmst)
      output_diff[l,]$rmst_diff_se <- sqrt(results[j,]$rmst_var + results[i,]$rmst_var)
      output_diff[l,]$rmst_diff_low <- output_diff[l,]$rmst_diff - qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_upp <- output_diff[l,]$rmst_diff + qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_pval <- 2*(1-pnorm(abs(output_diff[l,]$rmst_diff)/output_diff[l,]$rmst_diff_se))
      
      #--- RMST Ratio ---
      output_ratio[l,]$label_ratio <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_ratio[l,]$rmst_logratio <- (log(results[j,]$rmst) - log(results[i,]$rmst))
      output_ratio[l,]$rmst_logratio_se <- sqrt(results[j,]$rmst_var/(results[j,]$rmst^2) + results[i,]$rmst_var/(results[i,]$rmst^2))
      output_ratio[l,]$rmst_ratio <- exp(output_ratio[l,]$rmst_logratio)
      output_ratio[l,]$rmst_ratio_low <- exp(output_ratio[l,]$rmst_logratio - qnorm(1-alpha/2)*output_ratio[l,]$rmst_logratio_se)
      output_ratio[l,]$rmst_ratio_upp <- exp(output_ratio[l,]$rmst_logratio + qnorm(1-alpha/2)*output_ratio[l,]$rmst_logratio_se)
      output_ratio[l,]$rmst_logratio_pval <- 2*(1-pnorm(abs(output_ratio[l,]$rmst_logratio)/output_ratio[l,]$rmst_logratio_se))
      
      l <- l+1 #move to next row
    }
  }
  
  cat("\n\n\n")
  cat(paste('RMST calculated up to tau =',round(results$tau[1],3)))
  cat("\n\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) per Group \n\n")
  colnames(results) <- c("Group", "RMST", "Var", "SE", "Tau")
  rownames(results) <- c(paste("Group", results$Group,' '))
  print(round(results[c(2,4)],3))
  cat("\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) Differences \n\n")
  colnames(output_diff) <- c("Groups", "Est.", "SE", "CIL", "CIU", "p")
  rownames(output_diff) <- c(output_diff$Groups)
  print(round(output_diff[c(2,3,4,5,6)],3))
  cat("\n\n")
  
  cat ("Restricted Mean Survival Time (RMST) Ratios \n\n")  
  colnames(output_ratio) <- c("Groups", "Log Est.", "SE", "Est.", "CIL", "CIU", "p")
  rownames(output_ratio) <- c(output_ratio$Groups)
  print(round(output_ratio[c(2,3,4,5,6,7)],3))
}

