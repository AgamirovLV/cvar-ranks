
 let curdat_f,curdat_a,curdat_delta,curdat_betax,curdat_xu,curdat_n,curdat_r,curdat_ftype;

let time = performance.now();
 
//**********************процентные точки распределения коэффициента вариации************************************************
function datacvar(iflag) {

     c=window.open("","Cvar","toolbar=yes,menubar=yes,scroolbar=yes,width=650,height=400,left=400, top=150");
//Input:

//вектор вероятностей распределения
     let xw=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];

//вектор генеральных значений коэффициента вариации
     let gammcur=[0.05,0.3,0.5];

// минимальный объем выборки
     let nmin=3;

// максимальный объем выборки
     let nmax=10;

// верхний предел интегрирования
     curdat_xu=5.;

     let znaki=7;
     let nx=1;
     let stepx,lim,eps,n,gamm,tapp,cfvar;
     
 gammcur.forEach(function(gamm) {
     n=nmin-1;

 while(n<nmax) {
      n=n+1;curdat_f=n-1;
      c.document.write(n+" ");

//параметр нецентральности, если iflag=2
      curdat_delta=Math.sqrt(n)/gamm;

  xw.forEach(function(z) {
      curdat_betax=z;
      let x=[];
      
      if(iflag==0) {  
      x[1]=gamm;
      stepx=0.78;lim=1500;eps=0.000001;
      curdat_ftype="cvar";
//**********Nalder**********************************
      iter=simpl(x,nx,stepx,eps,lim,fun_cvar_t); //минимизация с интегрированием
//***************************************************
        //c=window.open("","Cvar","toolbar=yes,menubar=yes,scroolbar=yes,width=650,height=400,left=400, top=150");
        c.document.write(x[1].toFixed(znaki)+" ");
      }
     if(iflag==1) { //McKay’s approximation
        let zb=invchisquaredistribution(curdat_f,1-curdat_betax);
        gammat=Math.sqrt(zb/((1+1/(gamm*gamm))*(n-1)-zb*(n-1)/n))/gamm;
        c.document.write(gammat.toFixed(znaki)+" ");
      }
      if(iflag==2 || iflag==3 || iflag==4) { //распределение коэф.вариации через нецентральное t-распределение
         curdat_betax=1-curdat_betax;
         if(n<10) {
           tapp=curdat_delta+invnormaldistribution(curdat_betax);
        }
         else {
           tapp=invnontapp(curdat_betax,curdat_f,curdat_delta);
       }
          x[1]=tapp; //начальное приближение
          stepx=0.78;lim=1500;eps=0.000001;
          curdat_ftype="t";
//**********Nalder**********************************
         if(iflag==2) iter=simpl(x,nx,stepx,eps,lim,fun_cvar_t); //минимизация с интегрированием
         if(iflag==3) iter=simpl(x,nx,stepx,eps,lim,funtquantile_prncst); //минимизация с аппроксимацией ф.р. prncst
	 if(iflag==4) iter=simpl(x,nx,stepx,eps,lim,funtquantile_sf54r); //минимизация с аппроксимацией ф.р. sf54r
//***************************************************
         cfvar=curdat_delta/x[1];
         c.document.write(cfvar.toFixed(znaki)+" ");
      }
  });  //xw
    c.document.write("<br>");
 } //n
    c.document.write("<br>");c.document.write("<br>");

}); //gamm

   time = performance.now() - time;
   c.document.write("Program execution time=",(time/1000).toFixed(2)+"sec");

}
//**********************процентные точки нецентрального распределения Стьюдента*************************************************
function datatquantile(iflag) {

  c=window.open("","Student","toolbar=yes,menubar=yes,scroolbar=yes,width=650,height=400,left=400, top=150");
//Input:

//вектор вероятностей распределения
     let xw=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];

//вектор доверительных интервалов
     let xb=[0.99,0.95,0.9];

// минимальный объем выборки
     let nmin=3;

// максимальный объем выборки
     let nmax=10;

// верхний предел интегрирования
     curdat_xu=5.;

     let nx=1;
     let betap,zp,n;
     let stepx,lim,eps,tapp;
     let znaki=7;
 
   xb.forEach(function(z) {
       curdat_betax=z;
       n=nmin-1;
   
  while(n<nmax) { 
      n=n+1;curdat_f=n-1;
      c.document.write(n+" ");

    xw.forEach(function(betap) {
      
      stepx=0.18;lim=1500;eps=0.000001;
      curdat_ftype="t";
      let x=[];
      zp=invnormaldistribution(betap);
      curdat_delta=zp*Math.sqrt(n);  //параметр нецентральности может быть изменен пользователем
      if(n<10) {
       tapp=curdat_delta+invnormaldistribution(curdat_betax);
      }
      else {
       tapp=invnontapp(curdat_betax,curdat_f,curdat_delta);
      }
      x[1]=tapp;
//*******************Nalder*************************************
      if(iflag==0) iter=simpl(x,nx,stepx,eps,lim,fun_cvar_t); //минимизация с интегрированием
      if(iflag==1) iter=simpl(x,nx,stepx,eps,lim,funtquantile_prncst); //минимизация с аппроксимацией ф.р. prncst
      if(iflag==2) iter=simpl(x,nx,stepx,eps,lim,funtquantile_sf54r);  //минимизация с аппроксимацией ф.р. sf54r
//**************************************************************
       c.document.write(x[1].toFixed(znaki)+" ");
  }); //xw
    c.document.write("<br>");
 } //n
   c.document.write("<br>");c.document.write("<br>");
}); //xb
   time = performance.now() - time;
   c.document.write("Program execution time=",(time/1000).toFixed(2)+"sec");
}
//************************Числовые характеристики порядковых статистик (мат.ожидания и дисперсии)***********************
function dataorder(iflag) {
     c=window.open("","Order","toolbar=yes,menubar=yes,scroolbar=yes,width=650,height=400,left=400, top=150");
     let znaki=7;
     curdat_n=19; //объем выборки
     let alpha_nr=0,var_r,d;
     let r,rmax;
     let vorder=[];
     r=0; 
     if(iflag==0 || iflag==1) r=Math.floor(curdat_n/2); //для нормального распр. 
     rmax=curdat_n;
     n=curdat_n;
  while(r<rmax) { 
      r=r+1; 
      curdat_r=r;
      c.document.write(curdat_r+" ");
    if(iflag==0) { //прямое интегрирование, нормальное распределение
      curdat_ftype="Normal";
      curdat_order="mean";
      alpha_nr=IntegrateFunction("5",721,-9,9,forder);
      d=s2pi*gamma(r)*gamma(n-r+1)/gamma(n+1);
      alpha_nr=alpha_nr/d;
      c.document.write(alpha_nr.toFixed(znaki)+" ");
      curdat_order="var";
      var_r=IntegrateFunction("5",721,-9,9,forder);
      var_r=var_r/d-alpha_nr*alpha_nr;
      c.document.write(var_r.toFixed(znaki)+" ");
      c.document.write("<br>");
    }
   if(iflag==1) { //аппроксимация Дэйвида-Джонсона, нормальное распределение
      vorder=ordern(n,r/(n+1.),r/(n+1.));
      c.document.write(vorder[0].toFixed(znaki)+" ");
      c.document.write(vorder[1].toFixed(znaki)+" ");
      c.document.write("<br>");
    }
    if(iflag==2) {
      curdat_ftype="Weibull";
      curdat_order="mean";
      alpha_nr=IntegrateFunction("5",721,-15,15,forder); //прямое интегрирование, распределение Вейбулла
      d=gamma(r)*gamma(n-r+1)/gamma(n+1);
      alpha_nr=alpha_nr/d;
      c.document.write(alpha_nr.toFixed(znaki)+" ");
      curdat_order="var";
      var_r=IntegrateFunction("5",721,-15,15,forder);
      var_r=var_r/d-alpha_nr*alpha_nr;
      c.document.write(var_r.toFixed(znaki)+" ");
      c.document.write("<br>");
    }
 if(iflag==3) { //аппроксимация Дэйвида-Джонсона, распределение Вейбулла
      vorder=orderw(n,r/(n+1.),r/(n+1.));
      c.document.write(vorder[0].toFixed(znaki)+" ");
      c.document.write(vorder[1].toFixed(znaki)+" ");
      c.document.write("<br>");
    }
  }
}
//*******************интегрируемая функция t-распределения и коэф.вариации***************
function fpolinom_cvar_t(x) {
  let p,ds,z;
  a=curdat_a;
  delta=curdat_delta;
  f=curdat_f;
  ds=2*Math.pow(x,f-1)*Math.pow(f/2,f/2)*Math.exp(-f*x*x/2)/gamma(f/2); //плотность распределения с.к.о.
  if(curdat_ftype=="cvar") z=(x/a-1)*delta; 
  if(curdat_ftype=="t") z=x*a-delta;
  p =ds*normaldistribution(z);
  return(p);
}
//*****************Минимизируемые функции**********************************
function fun_cvar_t(x) {
    curdat_a=x[1];
    betax=curdat_betax;
    betar=IntegrateFunction("5",801,0,curdat_xu,fpolinom_cvar_t);
    if(curdat_ftype=="cvar") betar=1-betar; 
    return((betax-betar)*(betax-betar));
}
//********************************************************************
function funtquantile_prncst(x) {
    betax=curdat_betax;
    betar=prncst(x[1],curdat_f,curdat_delta);
    return((betax-betar)*(betax-betar));
}
//********************************************************************
function funtquantile_sf54r(x) {
    betax=curdat_betax;
    betar=sf54r(x[1],curdat_delta,curdat_f);
    return((betax-betar)*(betax-betar));
}
//**интегрируемая функция мат.ожидание и дисперсия порядковой статистики**
function forder(x) {
  let fr,r,n,cdf,df;
  r=curdat_r;
  n=curdat_n;
  if(curdat_ftype=="Normal") {
     cdf=normaldistribution(x);
     df=Math.exp(-x*x*0.5);
  }
  if(curdat_ftype=="Weibull") {
     cdf=1-Math.exp(-Math.exp(x));
     df=(1-cdf)*Math.exp(x);   
  }
     if(curdat_order=="mean") fr=x*df*(1-cdf)**(n-r)*cdf**(r-1);
     if(curdat_order=="var") fr=x*x*df*(1-cdf)**(n-r)*cdf**(r-1);
     return(fr);
}


