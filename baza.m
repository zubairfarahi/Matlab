function r = baza(sursa,b1,b2)
  if b1 < 2 || b1 > 30
    error('Base b1 need to be bigger than 1 and less than 31');
  else if b2 < 2 || b2 > 30
    error('Base b2 need to be bigger than 1 and less than 31');
  else
    bTranzition = base2dec(sursa,b1);
    r = dec2base(bTranzition,b2);
  end
end