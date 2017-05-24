function new_outputs=add_outputs(o1,o2)


new_outputs=o1;
fields=fieldnames(o1);
for ii=1:length(fields)
    f=new_outputs.(fields{ii});
    f2=o2.(fields{ii});
    if isstruct(f)
        f.fast=f.fast+f2.fast;
        f.slow=f.slow+f2.slow;
        f.deadmic=f.deadmic+f2.deadmic;
        new_outputs.(fields{ii})=f;
    else
        new_outputs.(fields{ii})=f+f2;
    end
end

end