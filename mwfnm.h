#ifdef _WIN32
# ifdef fnm_EXPORTS 
 /* We are building this library */
#  define FNM_EXPORT __declspec(dllexport)
# else
 /* We are using this library */
#  define FNM_EXPORT __declspec(dllimport)
# endif
#else
#  define FNM_EXPORT
#endif

FNM_EXPORT int testMe(int k);

FNM_EXPORT float useArray(const float* pData, const int nData);
