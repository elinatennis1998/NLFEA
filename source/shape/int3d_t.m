function [w,ss] = int3d_t(l,lint,ib)

% C-----------------------------------------------------------------------
%       subroutine int3d_t(ll,lint,s)
% 
% c      Purpose: Gauss quadrature for 3-d Tetrahedral element
% c      Inputs:
% c         ll       - Type of quadrature
% c      Outputs:
% c         lint     - Number of quadrature points
% c         s(4,*)   - Values of volume coordinates and weights
% c-----[--.----+----.----+----.-----------------------------------------]

% Modified 8/23/2011 for higher precision

   ss = zeros(3,1);
   s = zeros(4,16);
r9a = 0.0651301029022d00;
w9a = 0.0533472356088d00;
r9b = 0.8697397941956d00;
r9c = 0.3128654960049d00;
w9b = 0.0771137608903d00;
r9d = 0.6384441885698d00;
r9e = 0.0486903154253d00;
r9f = 0.2603459660790d00;
w9c = 0.1756152574332d00;
r9g = 0.4793080678419d00;
w9d = -0.1495700444677d00;
r9h = 0.3333333333333d00;
r7a = 0.33333333333333333333333333d00;
w7a = 0.225d00;
r7b = 0.059715871789769781d00;
w7b = 0.1323941527885061807376493878331519d00;
r7c = 0.4701420641051150897704412095134476d00;
r7d = 0.79742698535308731d00;
w7d = 0.1259391805448271525956839455001813d00;
r7e = 0.1012865073234563388009873619151238d00;
two = 2.d0;
r2 = 1/6;
w2 = 1/3;
   
   if ib == 0
       
% c     1 pt. quadrature O(h^2)
      if (lint == 1)

        for i = 1:3
          s(i,1) = 0.25d0;
        end
        s(4,1) = 1.0d0/6.d0;
        
      end
      
% c     4 pt. quadrature O(h^3)
      if(lint == 4) %verified 8/23/2011

        s(4,4) = 0.25d0/6.d0;
        for i = 1:3
          for j = 1:4
            s(i,j) = 0.1381966011250105151795413165634361;
          end
          s(i,i) = 0.58541019662496852;
          s(4,i) = s(4,4);
        end
      end
            
% c     5 pt quadrature
          if(lint == 5) %verified 8/23/2011

           s(1,1)=0.25d0;
           s(2,1)=0.25d0;
           s(3,1)=0.25d0;
           s(4,1)=-4.d0/30.d0;
           for i=2:5
           for j=1:3
             s(j,i)=1.d0/6.d0;
             if((i-1)==j)
             s(j,i)=0.5d0;
             end
           end
           s(4,i)=0.45d0/6.d0;
           end
          end

% c     11 pt. quadrature O(h^4)
      if(lint == 11) %verified 8/23/2011

%         for i = 1:3
%           for j = 1:10
%             s(i,j) = 0.0d0;
%           end
%           s(i, i  ) = 1.00d0;
%           s(i, i+4) = 0.50d0;
%           s(i, i+7) = 0.50d0;
%           s(i, 11 ) = 0.25d0;
%         end
%         s(2, 5) = 0.50d0;
%         s(3, 6) = 0.50d0;
%         s(1,10) = 0.50d0;
%         for j = 1:4
%           s(4,j) = 1.d0/360.d0;
%         end
%         for j = 5:10
%           s(4,j) = 1.d0/90.d0;
%         end
%         s(4,11) = 4.d0/45.d0;

          s = [0.25d0 0.25d0 0.25d0 -0.013155555555555555555555
               0.7857142857142857143d0 0.0714285714285714286d0 0.0714285714285714286d0 0.0076222222222222222222222222222222222d0
               0.0714285714285714286d0 0.7857142857142857143d0 0.0714285714285714286d0 0.0076222222222222222222222222222222222d0
               0.0714285714285714286d0 0.0714285714285714286d0 0.7857142857142857143d0 0.0076222222222222222222222222222222222d0
               0.0714285714285714286d0 0.0714285714285714286d0 0.0714285714285714286d0 0.0076222222222222222222222222222222222d0
               0.399403576166799d0 0.399403576166799d0 0.100596423833201d0 0.024888888888888888888888888888888888d0
               0.100596423833201d0 0.399403576166799d0 0.399403576166799d0 0.024888888888888888888888888888888888d0
               0.399403576166799d0 0.100596423833201d0 0.399403576166799d0 0.024888888888888888888888888888888888d0
               0.100596423833201d0 0.399403576166799d0 0.100596423833201d0 0.024888888888888888888888888888888888d0
               0.399403576166799d0 0.100596423833201d0 0.100596423833201d0 0.024888888888888888888888888888888888d0
               0.100596423833201d0 0.100596423833201d0 0.399403576166799d0 0.024888888888888888888888888888888888d0]';

      end

      if lint == 14 %verified 8/23/2011
          
         s = [0.3108859192633006097973457337634578 0.3108859192633006097973457337634578 0.3108859192633006097973457337634578 0.1126879257180158507991856523332863/6
              0.067342242210098213 0.3108859192633006097973457337634578 0.3108859192633006097973457337634578 0.1126879257180158507991856523332863/6
              0.3108859192633006097973457337634578 0.067342242210098213 0.3108859192633006097973457337634578 0.1126879257180158507991856523332863/6
              0.3108859192633006097973457337634578 0.3108859192633006097973457337634578 0.067342242210098213 0.1126879257180158507991856523332863/6
              
              0.0927352503108912264023239137370306 0.0927352503108912264023239137370306 0.0927352503108912264023239137370306 0.0734930431163619495437102054863275/6
              0.72179424906732637 0.0927352503108912264023239137370306 0.0927352503108912264023239137370306 0.0734930431163619495437102054863275/6
              0.0927352503108912264023239137370306 0.72179424906732637 0.0927352503108912264023239137370306 0.0734930431163619495437102054863275/6
              0.0927352503108912264023239137370306 0.0927352503108912264023239137370306 0.72179424906732637 0.0734930431163619495437102054863275/6
              
              0.0455037041256496494918805262793394 0.0455037041256496494918805262793394 0.45449629587435036 0.0425460207770814664380694281202574/6
              0.45449629587435036 0.0455037041256496494918805262793394 0.0455037041256496494918805262793394 0.0425460207770814664380694281202574/6
              0.0455037041256496494918805262793394 0.45449629587435036 0.0455037041256496494918805262793394 0.0425460207770814664380694281202574/6
              
              0.45449629587435036 0.0455037041256496494918805262793394 0.45449629587435036 0.0425460207770814664380694281202574/6
              0.0455037041256496494918805262793394 0.45449629587435036 0.45449629587435036 0.0425460207770814664380694281202574/6
              0.45449629587435036 0.45449629587435036 0.0455037041256496494918805262793394 0.0425460207770814664380694281202574/6]';
          
      end
      
% c     15 pt. quadrature O(h^5)
        if(lint == 15) %verified 8/23/2011
            
            s = [0.25d0 0.25d0 0.25d0 0.030283678097089d0
               0.0d0 0.333333333333333333333d0 0.333333333333333333333d0 0.006026785714286d0
               0.333333333333333333333d0 0.0d0 0.333333333333333333333d0 0.006026785714286d0
               0.333333333333333333333d0 0.333333333333333333333d0 0.0d0 0.006026785714286d0
               0.333333333333333333333d0 0.333333333333333333333d0 0.333333333333333333333d0 0.006026785714286d0
               0.727272727272727272727d0 0.090909090909090909090d0 0.090909090909090909090d0 0.011645249086029d0
               0.090909090909090909090d0 0.727272727272727272727d0 0.090909090909090909090d0 0.011645249086029d0
               0.090909090909090909090d0 0.090909090909090909090d0 0.727272727272727272727d0 0.011645249086029d0
               0.090909090909090909090d0 0.090909090909090909090d0 0.090909090909090909090d0 0.011645249086029d0
               0.066550153573664d0 0.066550153573664d0 0.433449846426336d0 0.010949141561386d0
               0.066550153573664d0 0.433449846426336d0 0.066550153573664d0 0.010949141561386d0
               0.433449846426336d0 0.066550153573664d0 0.066550153573664d0 0.010949141561386d0
               0.066550153573664d0 0.433449846426336d0 0.433449846426336d0 0.010949141561386d0
               0.433449846426336d0 0.066550153573664d0 0.433449846426336d0 0.010949141561386d0
               0.433449846426336d0 0.433449846426336d0 0.066550153573664d0 0.010949141561386d0]';
           
        end
        
% c     16 pt. quadrature O(h^5) % QUESTIONABLE accurate only to 8 places
        if(lint == 16)

        s(4,4) = 0.8395632516687135d-02;
        for i = 1:3
          for j = 1:4
            s(i,j) = 0.7611903264425430d-01;
          end
          s(i,i) = 0.7716429020672371d+00;
          s(4,i) = s(4,4);
        end
        for i = 5:16
          s(4,i) = 0.1109034477221540d-01;
        end

        s(1, 5) = 0.1197005277978019d+00;
        s(2, 5) = 0.7183164526766925d-01;
        s(3, 5) = 0.4042339134672644d+00;
        s(1, 6) = 0.4042339134672644d+00;
        s(2, 6) = 0.1197005277978019d+00;
        s(3, 6) = 0.7183164526766925d-01;
        s(1, 7) = 0.4042339134672644d+00;
        s(2, 7) = 0.4042339134672644d+00;
        s(3, 7) = 0.1197005277978019d+00;
        s(1, 8) = 0.7183164526766925d-01;
        s(2, 8) = 0.4042339134672644d+00;
        s(3, 8) = 0.4042339134672644d+00;

        s(1, 9) = 0.1197005277978019d+00;
        s(2, 9) = 0.4042339134672644d+00;
        s(3, 9) = 0.7183164526766925d-01;
        s(1,10) = 0.4042339134672644d+00;
        s(2,10) = 0.1197005277978019d+00;
        s(3,10) = 0.4042339134672644d+00;
        s(1,11) = 0.7183164526766925d-01;
        s(2,11) = 0.4042339134672644d+00;
        s(3,11) = 0.1197005277978019d+00;
        s(1,12) = 0.4042339134672644d+00;
        s(2,12) = 0.7183164526766925d-01;
        s(3,12) = 0.4042339134672644d+00;

        s(1,13) = 0.1197005277978019d+00;
        s(2,13) = 0.4042339134672644d+00;
        s(3,13) = 0.4042339134672644d+00;
        s(1,14) = 0.7183164526766925d-01;
        s(2,14) = 0.1197005277978019d+00;
        s(3,14) = 0.4042339134672644d+00;
        s(1,15) = 0.4042339134672644d+00;
        s(2,15) = 0.7183164526766925d-01;
        s(3,15) = 0.1197005277978019d+00;
        s(1,16) = 0.4042339134672644d+00;
        s(2,16) = 0.4042339134672644d+00;
        s(3,16) = 0.7183164526766925d-01;
        end
        
       if lint == 24 %verified 8/23/2011
         
         s = [0.2146028712591520292888392193862850 0.2146028712591520292888392193862850 0.2146028712591520292888392193862850 0.0399227502581674920996906275574800/6
              0.35619138622254387 0.2146028712591520292888392193862850 0.2146028712591520292888392193862850 0.0399227502581674920996906275574800/6
              0.2146028712591520292888392193862850 0.35619138622254387 0.2146028712591520292888392193862850 0.0399227502581674920996906275574800/6
              0.2146028712591520292888392193862850 0.2146028712591520292888392193862850 0.35619138622254387 0.0399227502581674920996906275574800/6
              
              0.0406739585346113531155794489564101 0.0406739585346113531155794489564101 0.0406739585346113531155794489564101 0.0100772110553206429480132374459369/6
              0.87797812439616596 0.0406739585346113531155794489564101 0.0406739585346113531155794489564101 0.0100772110553206429480132374459369/6
              0.0406739585346113531155794489564101 0.87797812439616596 0.0406739585346113531155794489564101 0.0100772110553206429480132374459369/6
              0.0406739585346113531155794489564101 0.0406739585346113531155794489564101 0.87797812439616596 0.0100772110553206429480132374459369/6
              
              0.3223378901422755103439944707624921 0.3223378901422755103439944707624921 0.3223378901422755103439944707624921 0.0553571815436547220951532778537260/6
              0.03298632957317360 0.3223378901422755103439944707624921 0.3223378901422755103439944707624921 0.0553571815436547220951532778537260/6
              0.3223378901422755103439944707624921 0.03298632957317360 0.3223378901422755103439944707624921 0.0553571815436547220951532778537260/6
              0.3223378901422755103439944707624921 0.3223378901422755103439944707624921 0.03298632957317360 0.0553571815436547220951532778537260/6
              
              0.0636610018750175252992355276057270 0.0636610018750175252992355276057270 0.6030056647916491413674311390609397 0.0482142857142857142857142857142857/6
              0.6030056647916491413674311390609397 0.0636610018750175252992355276057270 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6
              0.0636610018750175252992355276057270 0.6030056647916491413674311390609397 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6
              
              0.0636610018750175252992355276057270 0.0636610018750175252992355276057270 0.26967233145831571 0.0482142857142857142857142857142857/6
              0.0636610018750175252992355276057270 0.26967233145831571 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6
              0.26967233145831571 0.0636610018750175252992355276057270 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6
              
              0.0636610018750175252992355276057270 0.6030056647916491413674311390609397 0.26967233145831571 0.0482142857142857142857142857142857/6
              0.26967233145831571 0.0636610018750175252992355276057270 0.6030056647916491413674311390609397 0.0482142857142857142857142857142857/6
              0.6030056647916491413674311390609397 0.26967233145831571 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6
              
              0.0636610018750175252992355276057270 0.26967233145831571 0.6030056647916491413674311390609397 0.0482142857142857142857142857142857/6
              0.6030056647916491413674311390609397 0.0636610018750175252992355276057270 0.26967233145831571 0.0482142857142857142857142857142857/6
              0.26967233145831571 0.6030056647916491413674311390609397 0.0636610018750175252992355276057270 0.0482142857142857142857142857142857/6]';
           
       end
        
       if lint == 36 %verified 8/23/2011
         
         s = [0.0406107071929452723515318677627212 0.0406107071929452723515318677627212 0.0406107071929452723515318677627212 0.0061834158394585176827283275896323/6
              0.87816787842116417 0.0406107071929452723515318677627212 0.0406107071929452723515318677627212 0.0061834158394585176827283275896323/6
              0.0406107071929452723515318677627212 0.87816787842116417 0.0406107071929452723515318677627212 0.0061834158394585176827283275896323/6
              0.0406107071929452723515318677627212 0.0406107071929452723515318677627212 0.87816787842116417 0.0061834158394585176827283275896323/6
              
              0.1787522026964984761546314943983834 0.1787522026964984761546314943983834 0.1787522026964984761546314943983834 0.0785146502738723588282424885149114/6
              0.4637433919105046 0.1787522026964984761546314943983834 0.1787522026964984761546314943983834 0.0785146502738723588282424885149114/6
              0.1787522026964984761546314943983834 0.4637433919105046 0.1787522026964984761546314943983834 0.0785146502738723588282424885149114/6
              0.1787522026964984761546314943983834 0.1787522026964984761546314943983834 0.4637433919105046 0.0785146502738723588282424885149114/6
              
              0.3249495905373373335715573286644841 0.3249495905373373335715573286644841 0.3249495905373373335715573286644841 0.0447395776143224792777362432057442/6
              0.025151228387988001 0.3249495905373373335715573286644841 0.3249495905373373335715573286644841 0.0447395776143224792777362432057442/6
              0.3249495905373373335715573286644841 0.025151228387988001 0.3249495905373373335715573286644841 0.0447395776143224792777362432057442/6
              0.3249495905373373335715573286644841 0.3249495905373373335715573286644841 0.025151228387988001 0.0447395776143224792777362432057442/6
              
              0.1340777379721611918326213913378565 0.1340777379721611918326213913378565 0.7270125070093171000000000000000000 0.0121651445922912935604045907106762/6
              0.1340777379721611918326213913378565 0.7270125070093171000000000000000000 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              0.7270125070093171000000000000000000 0.1340777379721611918326213913378565 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              
              0.1340777379721611918326213913378565 0.1340777379721611918326213913378565 0.0048320170463606038 0.0121651445922912935604045907106762/6
              0.1340777379721611918326213913378565 0.0048320170463606038 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              0.0048320170463606038 0.1340777379721611918326213913378565 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              
              0.1340777379721611918326213913378565 0.7270125070093171000000000000000000 0.0048320170463606038 0.0121651445922912935604045907106762/6
              0.0048320170463606038 0.1340777379721611918326213913378565 0.7270125070093171000000000000000000 0.0121651445922912935604045907106762/6
              0.7270125070093171000000000000000000 0.0048320170463606038 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              
              0.1340777379721611918326213913378565 0.0048320170463606038 0.7270125070093171000000000000000000 0.0121651445922912935604045907106762/6
              0.7270125070093171000000000000000000 0.1340777379721611918326213913378565 0.0048320170463606038 0.0121651445922912935604045907106762/6
              0.0048320170463606038 0.7270125070093171000000000000000000 0.1340777379721611918326213913378565 0.0121651445922912935604045907106762/6
              
              0.0560275404597284769655799958528421 0.0560275404597284769655799958528421 0.3265740998664049580757011076659178 0.0280223074984909211766930561858945/6
              0.0560275404597284769655799958528421 0.3265740998664049580757011076659178 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6
              0.3265740998664049580757011076659178 0.0560275404597284769655799958528421 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6
              
              0.0560275404597284769655799958528421 0.0560275404597284769655799958528421 0.56137081921413812 0.0280223074984909211766930561858945/6
              0.0560275404597284769655799958528421 0.56137081921413812 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6
              0.56137081921413812 0.0560275404597284769655799958528421 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6
              
              0.0560275404597284769655799958528421 0.3265740998664049580757011076659178 0.56137081921413812 0.0280223074984909211766930561858945/6
              0.56137081921413812 0.0560275404597284769655799958528421 0.3265740998664049580757011076659178 0.0280223074984909211766930561858945/6
              0.3265740998664049580757011076659178 0.56137081921413812 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6
              
              0.0560275404597284769655799958528421 0.56137081921413812 0.3265740998664049580757011076659178 0.0280223074984909211766930561858945/6
              0.3265740998664049580757011076659178 0.0560275404597284769655799958528421 0.56137081921413812 0.0280223074984909211766930561858945/6
              0.56137081921413812 0.3265740998664049580757011076659178 0.0560275404597284769655799958528421 0.0280223074984909211766930561858945/6]';
           
       end

       ss(:) = s(1:3,l);
       w = s(4,l);
        
   else
      if (lint == 1)

        for i = 1:2
          s(i,1) = 1/3;
        end
        w = 1.0d0/two;
        
      end
       
       if lint == 3
           
           w = w2/two;
           s(1,1)=r2;
           s(2,1)=r2;
           s(3,1)=0;
           s(1,2)=r2;
           s(2,2)=1-2*r2;
           s(3,2)=0;
           s(1,3)=1-2*r2;
           s(2,3)=r2;
           s(3,3)=0;
           
       end
    if lint == 7 
        if l == 1  
            w= w7a/two;
        end
        if l>=2 && l<=4
            w=w7b/two;
        end
        if l>=5 && l<=7
            w=w7d/two;
        end
        s(1,1)=r7a;
        s(2,1)=r7a;
        s(3,1)=r7a;
        for i=2:4
            s(1,i)=r7c;
            s(2,i)=r7c;
            s(3,i)=r7c;
        end
        s(1,2)=r7b;
        s(2,3)=r7b;
        s(3,4)=r7b;
        for i=5:7
            s(1,i)=r7e;
            s(2,i)=r7e;
            s(3,i)=r7e;
        end 
        s(1,5)=r7d;
        s(2,6)=r7d;
        s(3,7)=r7d;
    end
       
       if lint == 13 
        if (l >= 1 && l <= 3)
          w=w9a/two;
        end
        if (l >= 4 && l <= 9)
          w=w9b/two;
        end
        if (l >= 10 && l <= 12)
          w=w9c/two;
        end
        if (l == 13)
          w=w9d/two;
        end

        s(1,1)=r9a;
        s(1,2)=r9b;
        s(1,3)=r9a;
        s(1,4)=r9c;
        s(1,5)=r9d;
        s(1,6)=r9e;
        s(1,7)=r9d;
        s(1,8)=r9c;
        s(1,9)=r9e;
        s(1,10)=r9f;
        s(1,11)=r9g;
        s(1,12)=r9f;
        s(1,13)=r9h;

        s(2,1)=s(1,1);
        s(2,2)=s(1,1);
        s(2,3)=s(1,2);
        s(2,4)=s(1,6);
        s(2,5)=s(1,4);
        s(2,6)=s(1,5);
        s(2,7)=s(1,6);
        s(2,8)=s(1,5);
        s(2,9)=s(1,4);
        s(2,10)=s(1,10);
        s(2,11)=s(1,10);
        s(2,12)=s(1,11);
        s(2,13)=s(1,13);
       end
    
       ss(:) = s(1:3,l);
       
   end
      
end
