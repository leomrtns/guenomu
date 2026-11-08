#!/bin/sh
# Missing optional Python is a skipped check, never a passing unexecuted test.
if test "x$1" = 'x:'; then
  echo 'Python 3.6 or later is required for the output regression check.'
  exit 77
fi
exec "$@"
