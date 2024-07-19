#ifndef _USER_H_
#define _USER_H_

// user defined variables


// example variables used in Wii camera testing - replace with your own
// variables
#ifdef USERHOOK_VARIABLES

#ifndef WII_CAMERA
#define WII_CAMERA 0
#endif

#if WII_CAMERA == 1
WiiCamera           ircam;
int                 WiiRange=0;
int                 WiiRotation=0;
int                 WiiDisplacementX=0;
int                 WiiDisplacementY=0;
#endif  // WII_CAMERA

#endif  // USERHOOK_VARIABLES

#endif  // _USER_H_

