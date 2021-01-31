import math
import numpy as np 


def RA_to_hours(RA: list) -> float:
    '''Return hours given a list of hours, minutes and days.'''
    hours, minutes, seconds = RA[0], RA[1], RA[2]
    return hours + minutes/60 + seconds/(60*60)


def number_days_till_midnight(RA_hours: float) -> int:
    '''Return the number of days when object found at RA (in hours) is at midnight.'''
    number_of_days = 0
    if RA_hours == 12:
        return 0
    else:
        midnighthours = 12
        while midnighthours < RA_hours:
            midnighthours += 3.945/60
            number_of_days += 1
        return number_of_days-1

if __name__ == "__main__":
    from datetime import datetime  
    from datetime import timedelta 

    vernal_equinox_date = datetime.strptime("2000-03-20", "%Y-%m-%d")

    #Alpha Centauri
    object_hours_alpha = RA_to_hours([14, 39, 37])
    number_of_days_alpha = number_days_till_midnight(object_hours_alpha)
    print('Number of Days (Alpha Centauri): ', number_of_days_alpha)
    print("Date on midnight (Alpha Centauri)", vernal_equinox_date+timedelta(days=number_of_days_alpha))

    #GN-z11
    object_hours_gn = RA_to_hours([12, 36, 25.46])
    number_of_days_gn = number_days_till_midnight(object_hours_gn)
    print('Number of Days (GN-z11): ',number_of_days_gn)
    print("Date on midnight (GN-z11)", vernal_equinox_date+timedelta(days=number_of_days_gn))

    #SagA*
    object_hours_sag = RA_to_hours([17, 45, 40.041])
    number_of_days_sag = number_days_till_midnight(object_hours_sag)
    print('Number of Days (SagA*): ', number_of_days_sag)
    print("Date on midnight (SagA*)", vernal_equinox_date+timedelta(days=number_of_days_sag))

    #3C273
    object_hours_3C = RA_to_hours([12, 29, 6.7])
    number_of_days_3C = number_days_till_midnight(object_hours_3C)
    print('Number of Days (3C273): ',number_of_days_3C)
    print("Date on midnight (3C273)", vernal_equinox_date+timedelta(days=number_of_days_3C))
