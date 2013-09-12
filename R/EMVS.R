EMVS<-function(Y,
               X,
               v0s,
               v1 = NULL,
               type,
               beta_init,
               sigma_init,
               epsilon,
               temperature = NULL,
               Z,
               mu,
               Sigma,
               p,
               a,
               b,
               a_v1,
               b_v1,
               v1_g){
  .Call( "EMVS",Y,X,v0s,v1,type,
        beta_init,sigma_init,epsilon,
        temperature,Z,
        mu,Sigma,p,a,b,a_v1,b_v1,v1_g,PACKAGE = "EMVSpackage" )
}

