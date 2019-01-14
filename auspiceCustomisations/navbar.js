import React from "@libraries/react";
import styled from '@libraries/styled-components';

const logoPNG = require("./logo.png");

const NavBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
  height: 100%;
`;

const Title = styled.span`
  padding: 0px;
  color: ${(props) => props.theme.color};
  font-size: 20px;
  font-weight: 400;
  letter-spacing: 0.5rem;
`;

const NavBar = ({narrativeTitle, sidebar}) => {
  if (!sidebar) return null;
  return (
    <NavBarContainer>
      <img alt="splashPage" style={{padding: "5px 5px"}} width="40px" src={logoPNG}/>
      <Title href="/">
        {"A/H3N2 antigenicity"}
      </Title>
      <div style={{minWidth: 10}}/>
    </NavBarContainer>
  );
};

export default NavBar;
