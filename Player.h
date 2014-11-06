#ifndef PLAYER_H
#define PLAYER_H

struct Player {
	Player() : isHuman(false), isBlack(false) {}
	Player(bool isHuman, bool isBlack) : isHuman(isHuman), isBlack(isBlack) {}
	bool isHuman;
	bool isBlack;
};

#endif
