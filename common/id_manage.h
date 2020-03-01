#pragma once
#include "common.h"

extern IdType g_Id = 0;
extern void ResetId()
{
	g_Id = 0;
}
extern IdType NewId()
{
	return g_Id++;
}