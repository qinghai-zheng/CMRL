function  w_s = soft(sigma,thld)
%%%%function  w_s = soft(sigma,thld) soft threholding 
          w_s = sign(sigma).*max(0,abs(sigma)-thld); 
end