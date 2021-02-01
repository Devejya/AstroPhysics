from astropy.io import ascii
import matplotlib.pyplot as plt


def dRAvsdDec(filename: str, columns) -> None:
    ''' Plot the deltaRA vs delta Dec given a file containing the position offset in
    milliarcseconds, of a star as a function of time'''
    data =ascii.read(filename, names=columns)
    fig,ax=plt.subplots()
    ax.scatter(data['RA'],data['Dec'])
    #plt.axis([-430,400,-550,280])
    plt.xlabel('$\Delta$RA (mas)')
    plt.ylabel('$\Delta$Dec (mas)')
    plt.show()
    
    return None

if __name__=="__main__":
    filename = 'phy270_ass1_parallax.txt'
    columns = ['Days', 'RA', 'Dec']
    dRAvsdDec(filename, columns)
    