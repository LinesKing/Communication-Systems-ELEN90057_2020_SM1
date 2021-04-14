function M = QAM_sub_channel2(EsTxN0linear,Hlinear,Pe,Mrange)
    EsTx = 1; %% average energy per transimiting symbol EsTx = (M-1)/3*Eg
    N0 = EsTx/EsTxN0linear;

    H = Hlinear;
    
    Eg = EsTx*3./(Mrange-1); %% Eg
    d = sqrt(2.*Eg); %% dmin
    
    PER = QAM128(EsTx,N0,H,Mrange(1));
    if PER < Pe
        M = Mrange(1);
    else
        PER = QAM64(EsTx,N0,H,Mrange(2));
        if PER < Pe
            M = Mrange(2);
        else
            PER = QAM32(EsTx,N0,H,Mrange(3));
            if PER < Pe
                M = Mrange(3);
            else
                PER = QAM16(EsTx,N0,H,Mrange(4));
                if PER < Pe
                M = Mrange(4);
                else
                    PER = QAM8(EsTx,N0,H,Mrange(5));
                    if PER < Pe
                    M = Mrange(5);
                    else
                        PER = QAM4(EsTx,N0,H,Mrange(6));
                        if PER < Pe
                        M = Mrange(6);
                        else
                            PER = QAM2(EsTx,N0,H,Mrange(7));
                            if PER < Pe
                            M = Mrange(7);
                            else
                                M = Mrange(8);
                            end
                        end
                    end
                end
            end
        end
    end

    
end

