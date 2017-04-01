function [elf]=external_el(elb,d,ua,va,vfluxf)
global dte2
    elf= elb-dte2.*((DXF(AXB(d).*ua)+DYF(AYB(d).*va))-vfluxf);
    elf = bcond1(elf);
end