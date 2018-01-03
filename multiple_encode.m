  function x = multiple_encode(str)
  
      % ENCODE  Translate text to dots and dashes.
      % encode('text')
   
      x = '';
      str = upper(str);
      for k = 1:length(str);
         ch = str(k);
         % A blank in the text is worth three in the code.
         if ch == ' '
            x = [x '   '];
         else
            x = [x morse_encode(ch) ' '];
         end
      end
   
   end % encode