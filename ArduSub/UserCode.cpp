#include "Sub.h"
#include "UserVariables.h"
#include "string.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wframe-larger-than=25000"

#define BUFFER_LEN      (46)

AP_InertialSensor &_ins = AP::ins();           // 获取IMU数据
Compass &_compass = AP::compass();   // 获取磁力计数据
AP_GPS &_gps = AP::gps();           // 获取GPS数据
AP_HAL::UARTDriver* _serial_nav;         // 获取GPS数据
Vector3f delta_angle, delta_velocity, omg, fsf, mag;
float dangle_dt, magnorm;
Location gps_loc;
uint64_t time_us;

unsigned char serial_buffer[128] = {0xEB, 0x90};

#ifdef USERHOOK_INIT
void Sub::userhook_init()
{
    // put your initialisation code here
    // this will be called once at start-up
    phikf_app = AP_PhiKF(TS);
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "PhiKF App Init...");

    // &_ins = AP::ins();           // 获取IMU数据
    // Compass &_compass = AP::compass();   // 获取磁力计数据
    // AP_GPS &_gps = AP::gps();           // 获取GPS数据
    _serial_nav = AP::serialmanager().get_serial_by_id(4);           // 获取GPS数据
    //_serial4->begin(460800);
    //GCS_SEND_TEXT(MAV_SEVERITY_INFO, "serial4 begin");
    // Vector3f delta_angle, delta_velocity, omg, fsf, mag;
    // float dangle_dt, magnorm;
    // Location gps_loc;
}
#endif

#ifdef USERHOOK_FASTLOOP
void Sub::userhook_FastLoop()
{
    phikf_app.t_start++;
    // put your 100Hz code here
    // if (phikf_app.t_start <= 10000) {       // 等待10s
    //     return;
    // }

    time_us = AP_HAL::micros64();

    // Adding PhiKF processing
    _ins.get_delta_angle(delta_angle, dangle_dt);
    omg.x = delta_angle.y / dangle_dt;
    omg.y = delta_angle.x / dangle_dt;
    omg.z = -delta_angle.z / dangle_dt;
    _ins.get_delta_velocity(delta_velocity, dangle_dt);
    fsf.x = delta_velocity.y / dangle_dt;
    fsf.y = delta_velocity.x / dangle_dt;
    fsf.z = -delta_velocity.z / dangle_dt;
    if (_compass.available()) {
        const Vector3f &field_Ga = _compass.get_field();
        magnorm = sqrt(field_Ga.x*field_Ga.x + field_Ga.y*field_Ga.y + field_Ga.z*field_Ga.z);
        if(magnorm>1e-6) mag = field_Ga/magnorm;
    }
    
    phikf_app.setIMU(omg.x, omg.y, omg.z, fsf.x, fsf.y, fsf.z);
    phikf_app.setMag(mag.x, mag.y, mag.z);

    // GCS_SEND_TEXT(MAV_SEVERITY_INFO, "imu:%f %f %f %f %f %f", phikf_app.imub.wm.i, phikf_app.imub.wm.j, phikf_app.imub.wm.k, phikf_app.imub.vm.i, phikf_app.imub.vm.j, phikf_app.imub.vm.k);

    if (_gps.get_hdop()<=5000) {
        gps_loc = _gps.location();
        phikf_app.setGNSS(gps_loc.lat, gps_loc.lng, gps_loc.alt);
        //GCS_SEND_TEXT(MAV_SEVERITY_INFO, "g: %f %f %f", phikf_app.gps.pos.i, phikf_app.gps.pos.j, phikf_app.gps.pos.k);
    }
    
    memcpy(serial_buffer+2, &time_us, 8);
    // 不让用下边这种写法， [-Werror=cast-align] 垃圾警告
    // *(float*)(serial_buffer+10) = omg.y;  *(float*)(serial_buffer+14) = omg.x;  *(float*)(serial_buffer+18) = -omg.z;
    // *(float*)(serial_buffer+22) = fsf.y;  *(float*)(serial_buffer+26) = fsf.x;  *(float*)(serial_buffer+30) = -fsf.z;
    // *(float*)(serial_buffer+34) = mag.y;  *(float*)(serial_buffer+38) = mag.x;  *(float*)(serial_buffer+42) = -mag.z;
    memcpy(serial_buffer+10, &omg.x, 12);
    memcpy(serial_buffer+22, &fsf.x, 12);
    memcpy(serial_buffer+34, &mag.x, 12);

    _serial_nav->write(serial_buffer, BUFFER_LEN);

    #ifndef TURNOFF_NAV
    phikf_app.TimeUpdate();
    #endif
    // AP::logger().Write("IMUD", "TimeUS,state,wx,wy,wz,fx,fy,fz, mx,my,mz, lat,lng, alt", "QBffffffffffff",
    //                     time_us, phikf_app.state, phikf_app.imub.wm.i, phikf_app.imub.wm.j, phikf_app.imub.wm.k,
    //                     phikf_app.imub.vm.i, phikf_app.imub.vm.j, phikf_app.imub.vm.k,
    //                     phikf_app.magb.mag.i, phikf_app.magb.mag.j, phikf_app.magb.mag.k,
    //                     gps_loc.lat, gps_loc.lng, gps_loc.alt);
    // AP::logger().Write("PNAV", "TimeUS, pitch,roll,yaw, vE,vN,vU, Lat,Lng,Alt", "Qfffffffff",
    //                     time_us, phikf_app.ins_avp.att.i, phikf_app.ins_avp.att.j, phikf_app.ins_avp.att.k,
    //                     phikf_app.ins_avp.vn.i, phikf_app.ins_avp.vn.j, phikf_app.ins_avp.vn.k,
    //                     phikf_app.ins_avp.pos.i, phikf_app.ins_avp.pos.j, phikf_app.ins_avp.pos.k);

    //GCS_SEND_TEXT(MAV_SEVERITY_INFO, "Sreial4 send");
}
#endif

#ifdef USERHOOK_50HZLOOP
void Sub::userhook_50Hz()
{
    // put your 50Hz code here
}
#endif

#ifdef USERHOOK_MEDIUMLOOP
void Sub::userhook_MediumLoop()
{
    // put your 10Hz code here
}
#endif

#ifdef USERHOOK_SLOWLOOP
void Sub::userhook_SlowLoop()
{
    // put your 3.3Hz code here
}
#endif

#ifdef USERHOOK_SUPERSLOWLOOP
void Sub::userhook_SuperSlowLoop()
{
    // put your 1Hz code here
}
#endif

#pragma GCC diagnostic pop
