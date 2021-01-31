import math


def visibility_hours(A: float, delta: list, observatory: str, dark_hours: float) -> float:
    ''' '''
    if observatory == "GN":
        beta = 19.8238        
    if observatory == "GS":
        beta = -30.24075
    delta_radians = math.radians(delta[0] + delta[1]/60 + delta[2]/(60*60))
    beta_radians = math.radians(beta)
    HA = math.degrees(math.acos((math.sin(A)/(math.cos(delta_radians)*math.cos(beta_radians))) - 
        (math.tan(delta_radians)*math.tan(beta_radians))))

    return min(2*HA/6, dark_hours)


if __name__ == "__main__":
    A = math.radians(30)
    #Alpha Centauri
    alpha_dark_hours = 24-11.79
    print("Alpha Centauri Visibility Hours: ", visibility_hours(A, [-60, 50, 2], "GS", alpha_dark_hours))

    #GN-z11
    gn_dark_hours = 24-13.04
    print("GN-z11 Visibility Hours: ", visibility_hours(A, [62, 14, 31.4], "GN", gn_dark_hours))

    #SagA*
    sag_dark_hours = 24-11.08
    print("Sag-A* Visibility Hours: ", visibility_hours(A, [-29, 00, 28.1], "GS", sag_dark_hours))

    #3C273
    c_dark_hours = 24-13.00
    print("3C273 Visibility Hours: ", visibility_hours(A, [+2, 3, 9], "GN", c_dark_hours))