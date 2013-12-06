MatLab
======
function f=akin(n,x)
% Bickley function of order n >= 1 at x >= 0
% function f=akin(n,x)
%
   if n < 1 
      error('n must be > 0')
   elseif x < -1.0e-10
      error(sprintf('x must be >= 0. x=%.5g.',x))
   elseif x <= 0
      bick0=[ 1.57079632679489, 1.0, .78539816339745, .66666666666667, .58904862254808, .53333333333333, ...
      .4908738521234, .45714285714286, .42951462060798, .4063492063492 ] ;
      f=bick0(n) ;
   elseif x < 1
      xsq = x*x ;
      lnx = log(x) ;
      xi1 = x - .91169076342161 ;
      xi2 = x - .2451000192866 ;
      xi3 = x - 1.0 ;
      xi4 = x - .68448452306295 ;
      if n == 1
         f = (((((((((.7945473334662959e-4 * x + .51674609093834e-4) * x + .100122044600498e-1) * x + .8766262885587739e-2) ...
         * x + .5720051721276178) * x + .5454840912170553) * x + .139868293763185e2) * x + .1531257133183402e2) * x ...
         + .1501664584981108e3) * xi3 -.3916967515498982e2) / ((x + .1219596245756669e1) * xi3 -.1193155300614385e3) ...
         + (x * lnx * ((((.2337898258663651e-2 * xsq + .4646979722852471) * xsq + .3695696751512241e2) * xsq ...
         + .123465484355545e4) * xsq + .175237360009281e5) / ((xsq -.2256564898552151e3) * xsq + .175237360009281e5)) ;
      elseif n==2
         f = ((((.1403059e-10 * xi1 + .4811961706232723e3) + ((((.9084569646859357 * x + .490556407762756e2) * x ...
         + .7521727532834893e2) * x + .1693807200769639e4) + (((.1086764750096697e-1 * xi2 -1.0e-15) * xi2 ...
         + .1190096804348251e1) * (x^4))) * (xi1 * xi1)) / (.1810502008060146e4 + (((xi3 -.1893578319929816e2) * x ...
         -.7855291318496802e2) * xi3))) + (xsq * lnx * ((-.1520774316867189e9) / (((((xsq -.5631701819761997e2) * xsq ...
         -.3143123637802091e3) * xsq + .211218654889524e6) * xsq -.1267311930720872e8) * xsq + .3041548633734379e9)))) ;
      elseif n==3
         f = (((((.1173330390873767e5 * x + .1095667013274141e7) * xi4 -.379051e-8) * xi4 + .1269499275481224e7) ...
         * xi3 -.7043581454636306e6) /  (((((((xi3 -.5812262590904993e1) * x -.17551758398419e2) * x -.133516191424771e3) ...
         * x -.3942586515380026e4) * xi3 + .6077621585261822e5) * x + .2053835980116203e6) * xi3 -.2961415636470914e7)) ...
         + (x^3 * lnx * ((((.2631126488553487e-2 * xsq + .5562992588150486) * xsq + .3721363059831219e2) * xsq ...
         + .3027965765686327e4) / ((xsq -.2309130812632629e3) * xsq + .1816779459411791e5))) ;
      else
         ak1 = akin(n-1,x) ;
         ak2 = akin(n-2,x) ;
         ak3 = akin(n-3,x) ;
         f = (ak3-ak1)*x/(n-1)+(n-2)*ak2/(n-1) ;
      end
   elseif x < 6
      sqrtx = sqrt(x) ;
      expx = exp(-x) ;
      xrec = 1.0/x ;
      if n == 1
         f = (((((((.1822929159877549e2 * xrec + .3272001530672078e3) * xrec + .1326511766009986e4) * xrec ...
         + .1868734192859498e4) * xrec + .1059016416894119e4) * xrec + .2427580524508585e3) * xrec + .1833164538368226e2) ...
         * expx) / ((((((((xrec + .6590511376539962e2) *xrec + .5952592332227032e3) * xrec + .1687760486772990e4) * xrec ...
         + .1922624187690926e4) * xrec + .957005687628236e3) * xrec + .2028344160151355e3) * xrec + .1462653804563246e2) ...
         * sqrtx) ;
      elseif n==2
         f = ((((((((.8407469297501269e-1 * xrec + .5596498537189973e1) * xrec + .4801733781249936e2) * xrec ...
         + .123974074193467e3) * xrec + .12440411402683e3) * xrec + .5320941946830476e2) * xrec + .9534267279889207e1) ...
         * xrec + .5766817227841408) * expx) / ((((((((xrec + .2014290370371339e2) * xrec + .9747773947009136e2) * xrec ...
         + .1754089481769652e3) * xrec + .1383006201574071e3) * xrec + .5035544525458363e2) * xrec + .8124881079392082e1) ...
         * xrec + .4601255143693006) * sqrtx) ;
      elseif n==3
         f = ((((((((.6753131438170179e-1 * xrec + .4296855364707056e1) * xrec + .3523089103388543e2) * xrec ...
         + .8663638447767516e2) * xrec + .822191552088093e2) * xrec + .3295675793097689e2) * xrec + .5484114159179143e1) ...
         * xrec + .3057730071278506) * expx) / ((((((((xrec + .1898628899818902e2) * xrec + .8502251629078864e2) * xrec ...
         + .1397165642636965e3) * xrec + .9980505157913694e2) * xrec + .3279841354751966e2) * xrec + .4772136519722949e1) ...
         * xrec + .2439716682658748) * sqrtx) ;
      else
         ak1 = akin(n-1,x) ;
         ak2 = akin(n-2,x) ;
         ak3 = akin(n-3,x) ;
         f = (ak3-ak1)*x/(n-1)+(n-2)*ak2/(n-1) ;
      end
   elseif x < 673.5
      sqrtx = sqrt(x) ;
      expx = exp(-x) ;
      xrec = 1.0/x ;
      if n == 10
         f = (((((((((.2115457456609434e1 * xrec + .8380471295924849e1) * xrec + .8279846050989022e1) * xrec ...
         + .3176720369338028e1) * xrec + .5612690964007783) * xrec + .4842493199976428e-1) * xrec + .1966652996666221e-2) ...
         * xrec + .298485739150359e-4) / ((((((((xrec + .1303324507927971e2) * xrec + .2421477559533582e2) ...
         * xrec + .1531699714433142e2) * xrec + .4348767092235221e1) * xrec + .622422739384902) * xrec ...
         + .4627626853651856e-1 ) * xrec + .1691217608476688e-2) * xrec + .2381571628879027e-4) ) / sqrtx) * expx ;
      elseif n==9
         f = (((((((((.2271732886275686e1 * xrec + .9312166436472808e1) * xrec + .9546586415461775e1) * xrec ...
         + .3808588292608732e1) * xrec + .7008898871037985) * xrec + .6307055224810387e-1) * xrec + .2673933491122202e-2) ...
         * xrec + .4238334538408917e-4) / ((((((((xrec + .134693401036622e2) * xrec + .2599209729530514e2) * xrec ...
         + .1715435708724829e2) * xrec + .5098663790807535e1) * xrec + .7653242147414652) * xrec + .5971301111405558e-1) ...
         * xrec + .2289893952421851e-2) * xrec + .3381701691714022e-4) ) / sqrtx) * expx ;
      elseif n==8
         f = (((((((((.8986422681366424e-1 * xrec + .1112575911492456e1) * xrec + .1911501312353472e1) * xrec ...
         + .1072628589602992e1) * xrec + .2548459851246907) * xrec + .2815600768166693e-1) * xrec + .1418798411613728e-2) ...
         * xrec + .2613746296365794e-4) /  (((((((xrec + .384634427957306e1) * xrec + .3799408842694803e1) * xrec ...
         + .1507878566796656e1) * xrec + .2845497743872052) * xrec + .2689392859262307e-1) * xrec + .1218062894916834e-2) ...
         * xrec + .2085467815725935e-4)) / sqrtx) * expx ;
      else
         ak1 = akin(n+1,x) ;
         ak2 = akin(n+2,x) ;
         ak3 = akin(n+3,x) ;
         f = (((n+2)*ak3)-((n+1)*ak1))/x+ak2 ;
      end
   else
      f = 0.0 ;
   end
