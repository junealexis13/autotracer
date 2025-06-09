from datetime import datetime, timezone


def convert_unix_to_dt(unix_ts):
    """
    Convert a Unix timestamp in milliseconds to a UTC datetime object.
    
    :param unix_ts: Unix timestamp in milliseconds
    :return: UTC datetime object
    """
    # convert ms to s
    timestamp_s = unix_ts / 1000

    # conv to datetime (UTC)
    dt = datetime.fromtimestamp(timestamp_s, tz=timezone.utc)
    return dt
