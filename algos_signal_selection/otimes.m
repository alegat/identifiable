% a,b,c must belong to {-1, 0, 1}

function c = otimes(a,b)
    switch a
        case 1
            switch b
                case 1
                    c = 1;
                case 0
                    c = 0;
                case -1
                    c = -1; % only difference from odot
            end
        case 0
            switch b
                case 1
                    c = 0;
                case 0
                    c = 0;
                case -1
                    c = 0;
            end
        case -1
            switch b
                case 1
                    c = 1;
                case 0
                    c = 0;
                case -1
                    c = -1;
            end
    end
end