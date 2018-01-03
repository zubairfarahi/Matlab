function x = morse_encode(c)
      % ENCODE_CH  Translate one character to dots and dashes.
   
      S = {morse};
      X = {''};
      while ~isempty(S)
         N = S{1};
         x = X{1};
         S = S(2:end);
         X = X(2:end);
         if ~isempty(N)
            if N{1} == c;
               return
            else
               S = {N{2} N{3} S{:}};
               X = {[x '.'] [x '-'] X{:}};
            end
         end
      end
      x = '*';
   
   end % encode_ch
