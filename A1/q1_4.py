import math

def max_altitude(delta: list, observatory: str) -> float:
    ''' '''
    if observatory == "GN":
        beta = 19.8238
    if observatory == "GS":
        beta = -30.24075
    delta_radians = math.radians(delta[0] + delta[1]/60 + delta[2]/(60*60))
    beta_radians = math.radians(beta)
    return math.degrees(math.asin(math.sin(delta_radians)*math.sin(beta_radians) + math.cos(beta_radians)*math.cos(delta_radians)))


if __name__ == "__main__":
    #Alpha Centauri
    print ("Max Altitude of Alpha Centauri at GS: ", max_altitude([-60, 50, 2], "GS"))
    #GN-z11
    print ("Max Altitude of GN-z11 at GN: ", max_altitude([62, 14, 31.4], "GN"))
    #SagA*
    print ("Max Altitude of SagA* at GS: ", max_altitude([-29, 00, 28.1], "GS"))
    #3C273
    print ("Max Altitude of 3C273 at GN: ", max_altitude([+2, 3, 9], "GN"))


