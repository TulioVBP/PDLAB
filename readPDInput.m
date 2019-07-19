function readPDInput(input)
    fileID = fopen(input);
    tline = fgetl(fileID);
    while tline ~= "!END"
        b_sectionChange = strfind(tline,'#');
        if b_section
           switch tline
               case '# MATERIAL DATA'
               case '# PD DATA'
               case '# MESH DETAILS'
               case '# SIMULATION DETAILS'
               otherwise
                   error('Input file has an unknown section.')
           end
        end
        
    end
    
end