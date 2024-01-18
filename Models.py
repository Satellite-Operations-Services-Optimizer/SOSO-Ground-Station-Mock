from pydantic import BaseModel, Field
from datetime import datetime, timedelta
from typing import Optional

class schedule(BaseModel):
    id: int
    satellite_id:int
    ground_station_id: int
    asset_type: int
    start_time: datetime
    end_time: datetime
    status: str

       
class ground_station_schedule(BaseModel):
    id: int
    satellite_id:int
    ground_station_id: int
    asset_type: int
    start_time: datetime
    end_time: datetime
    status: str