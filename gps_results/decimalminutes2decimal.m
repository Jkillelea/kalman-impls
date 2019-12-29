function ddd = decimalminutes2decimal(decimalminutes)
    degrees = sign(decimalminutes) .* floor(abs(decimalminutes)/100);
    minutes = (decimalminutes - 100*degrees);
    decimal = minutes ./ 60;
    ddd = degrees + decimal;
end
